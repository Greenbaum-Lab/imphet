import os
import argparse
import pandas as pd
import numpy as np
import cupy as cp
import cupyx.scipy.special as cpx_spl
from bed_reader import open_bed


def init_bed_reader(bed_prefix, samples, threads):
	'''
	Return (bed_file, sample_indices, bim)

	sample_indices follow the order of `samples`.
	bim holds the chromosome / position columns loaded with pandas
	to avoid the slow np.loadtxt call inside bed_reader.
	'''
	bed_file = open_bed(f'{bed_prefix}.bed', num_threads = threads)
	bim = pd.read_csv(
		f'{bed_prefix}.bim',
		sep   = r'\s+',
		usecols = [0, 3],
		names = ['chr', 'pos'],
		dtype = {'chr': str, 'pos': np.int32}
	)
	iid = bed_file.iid
	common, idx_samples, idx_iid = np.intersect1d(samples, iid, return_indices = True)
	if len(common) < len(samples):
		raise KeyError('Some samples not found in FAM file')
	order = np.argsort(idx_samples)
	sample_indices = idx_iid[order].astype(np.int32)
	return bed_file, sample_indices, bim


def load_genotypes(bed_file, sample_indices, bim, region):
	chr = bim['chr'].values
	pos = bim['pos'].values
	variant_indices = np.where((chr == str(region[0])) & (pos >= region[1]) & (pos < region[2]))[0]
	genotypes = np.ascontiguousarray(bed_file.read(np.s_[sample_indices, variant_indices], dtype='int8'))
	return genotypes, pos[variant_indices]


def parse_labelled_regions(labelled_regions_str):
	labelled_regions = {}
	for region in labelled_regions_str.split(';'):
		label, region_str = region.split('=')
		chr_label, pos_range = region_str.split(':')
		start, end = map(int, pos_range.split('-'))
		labelled_regions[label] = [chr_label, start, end]
	return labelled_regions


def he_to_pe(heterozygosity):
	if heterozygosity > 0.5:
		return np.nan
	discriminant = 1 - 2 * heterozygosity
	return np.nan if discriminant < 0 else (1 - np.sqrt(discriminant)) / 2


def ha_update(allele_freq, selection_coeff):
	numerator = allele_freq * (1 - selection_coeff * allele_freq)
	denominator = 1 - selection_coeff * (allele_freq ** 2 + (1 - allele_freq) ** 2)
	return np.nan if abs(denominator) < 1e-12 else numerator / denominator


def iterate_ha_model(start_freq, selection_coeff, generations):
	freq = start_freq
	for _ in range(generations):
		freq = ha_update(freq, selection_coeff)
		if not (0 < freq < 1):
			return np.nan
	return freq


def bracketed_bisection(start_freq, end_freq, generations,
		lower_bound = -0.2, upper_bound = 0.2,
		tolerance = 1e-8, max_iterations = 1000):
	if abs(end_freq - start_freq) < tolerance:
		return 0.0

	def selection_difference(selection_coeff):
		return iterate_ha_model(start_freq, selection_coeff, generations) - end_freq

	left, right = lower_bound, upper_bound
	left_val, right_val = selection_difference(left), selection_difference(right)
	if np.isnan(left_val) or np.isnan(right_val) or left_val * right_val > 0:
		return np.nan

	for _ in range(max_iterations):
		mid = 0.5 * (left + right)
		mid_val = selection_difference(mid)
		if np.isnan(mid_val) or abs(mid_val) < tolerance or abs(right - left) < tolerance:
			return mid
		if left_val * mid_val > 0:
			left, left_val = mid, mid_val
		else:
			right, right_val = mid, mid_val
	return mid


def selection_coeff_null(het_counts_df, total_counts_df, sample_dates,
		time_window_start = 12000,
		time_window       = 2000,
		time_step         = 1000,
		generation_time   = 25):
	meta_cols = 3 if 'chr' in het_counts_df.columns else 1

	het_mat = het_counts_df.iloc[:, meta_cols:].to_numpy(dtype = np.float32)
	tot_mat = total_counts_df.iloc[:, meta_cols:].to_numpy(dtype = np.float32)
	with np.errstate(divide = 'ignore', invalid = 'ignore'):
		het_frac = np.where(tot_mat > 0, het_mat / tot_mat, np.nan)

	n_windows, _ = het_frac.shape
	temporal_starts = np.arange(0, time_window_start, time_step, dtype = int)[::-1]
	temporal_ends   = temporal_starts + time_window
	T_gen = int(time_step / generation_time)

	records = []
	for t in range(len(temporal_starts) - 1):
		s0, e0 = temporal_starts[t],   temporal_ends[t]
		s1, e1 = temporal_starts[t+1], temporal_ends[t+1]

		idx0 = np.where((sample_dates >= s0) & (sample_dates < e0))[0]
		idx1 = np.where((sample_dates >= s1) & (sample_dates < e1))[0]
		if idx0.size == 0 or idx1.size == 0:
			continue

		p0 = np.nanmean(het_frac[:, idx0], axis = 1)
		p1 = np.nanmean(het_frac[:, idx1], axis = 1)

		for w in range(n_windows):
			q0 = he_to_pe(p0[w])
			q1 = he_to_pe(p1[w])
			if np.isnan(q0) or np.isnan(q1):
				continue
			s_est = bracketed_bisection(q0, q1, T_gen)
			if np.isnan(s_est):
				continue
			s_est = max(0.0, s_est)
			meta_vals = het_counts_df.iloc[w, :meta_cols].tolist()
			records.append(meta_vals + [s0, e1, s_est])

	if 'chr' in het_counts_df.columns:
		col_hdr = ['chr', 'start', 'end', 'time_start', 'time_end', 's_coeff']
	else:
		col_hdr = ['label', 'time_start', 'time_end', 's_coeff']
	return pd.DataFrame.from_records(records, columns = col_hdr)


def compute_het(bed_file, sample_bed_indices, bim, samples, regions,
		batch_size = 100000, snp_density = 30):
	'''
	Compute per-region heterozygous and total SNP counts.
	'''
	n_regions  = len(regions)
	n_samples  = len(samples)
	max_size   = batch_size * snp_density

	het_counts        = np.empty((n_regions, n_samples), dtype = np.int32)
	total_snp_counts  = np.empty_like(het_counts)

	idx = 0
	while idx < n_regions:
		chr   = regions.at[idx, 'chr']
		start_pos = regions.at[idx, 'start']

		in_chr    = np.where(regions['chr'] == chr)[0]
		in_block  = in_chr[in_chr >= idx]
		block_end_candidates = in_block[regions.loc[in_block, 'end'] <= start_pos + max_size]
		end_region = block_end_candidates[-1]
		end_pos    = regions.at[end_region, 'end']

		genotypes, positions = load_genotypes(
			bed_file, sample_bed_indices, bim,
			[chr, int(start_pos), int(end_pos)]
		)

		for region_idx in range(idx, end_region + 1):
			region = regions.loc[region_idx]
			var_idx = np.where(
				(positions >= region['start']) & (positions < region['end'])
			)[0]
			g = genotypes[:, var_idx]
			het_counts[region_idx]       = (g == 1   ).sum(axis = 1)
			total_snp_counts[region_idx] = (g != -127).sum(axis = 1)

		idx = end_region + 1

	het_df = pd.concat(
		[regions.reset_index(drop = True),
		 pd.DataFrame(het_counts, columns = samples)],
		axis = 1
	)
	total_snp_df = pd.concat(
		[regions.reset_index(drop = True),
		 pd.DataFrame(total_snp_counts, columns = samples)],
		axis = 1
	)
	return het_df, total_snp_df


def compute_window_het(bed_file, sample_bed_indices, bim, samples,
		window_size = 1_000_000,
		step_size = 1_000_000,
		exclude_regions = None,
		block_size = 100_000):
	'''
	Return two DataFrames (het_counts_df, total_counts_df) that hold per-sample
	SNP counts for fixed genomic windows. Both share the columns chr, start, end, <one column per sample>
	'''
	col_hdr = ['chr', 'start', 'end'] + list(samples)
	het_rows = []
	total_rows = []

	bim_chr = bim['chr'].values
	bim_pos = bim['pos'].values

	for chr in pd.unique(bim_chr):
		pos_full = bim_pos
		chr_mask = bim_chr == chr
		chr_pos = pos_full[chr_mask]
		if chr_pos.size == 0:
			continue

		win_starts = np.arange(0, chr_pos.max() + 1, step_size, dtype = int)
		het_mat = np.zeros((win_starts.size, len(samples)), dtype = np.int32)
		cov_mat = np.zeros_like(het_mat)

		var_idx = np.where(chr_mask)[0]
		if exclude_regions:
			exclude = np.zeros(var_idx.size, dtype = bool)
			for ex_chr, ex_s, ex_e in exclude_regions:
				if ex_chr != chr:
					continue
				exclude |= (chr_pos >= ex_s) & (chr_pos < ex_e)
			var_idx = var_idx[~exclude]

		for b in range(0, var_idx.size, block_size):
			block = var_idx[b : b + block_size]
			g = bed_file.read(np.s_[sample_bed_indices, block], dtype = 'int8')
			win_idx = (pos_full[block] // step_size).astype(int)
			for k, w in enumerate(win_idx):
				het_mat[w] += (g[:, k] == 1)
				cov_mat[w] += (g[:, k] != -127)

		for i, ws in enumerate(win_starts):
			het_rows.append([chr, ws, ws + window_size] + het_mat[i].tolist())
			total_rows.append([chr, ws, ws + window_size] + cov_mat[i].tolist())

	het_df = pd.DataFrame(het_rows, columns = col_hdr)
	total_df = pd.DataFrame(total_rows, columns = col_hdr)
	return het_df, total_df


def compute_labelled_counts(bed_file, sample_bed_indices, bim, samples, labelled_regions):
	labels = list(labelled_regions.keys())
	n_samples = len(samples)
	het_arr = np.zeros((len(labels), n_samples), dtype = np.int32)
	total_arr = np.zeros_like(het_arr)

	for i, label in enumerate(labels):
		chr_label, start, end = labelled_regions[label]
		genotypes, _ = load_genotypes(bed_file, sample_bed_indices, bim, [chr_label, start, end])
		het_arr[i] = (genotypes == 1).sum(axis = 1)
		total_arr[i] = (genotypes != -127).sum(axis = 1)

	het_df = pd.DataFrame(het_arr,   columns = samples)
	total_df = pd.DataFrame(total_arr, columns = samples)
	het_df.insert(0, 'label', labels)
	total_df.insert(0, 'label', labels)
	return het_df, total_df


def compute_labelled_het(het_counts_df, total_counts_df, samples):
	frac = het_counts_df.iloc[:, 1:].div(total_counts_df.iloc[:, 1:])
	labels = het_counts_df['label'].to_numpy()
	frac = frac.T
	frac.columns = labels
	frac.index = samples
	return frac


def cupy_parallel_linregress(x, y):
	'''
	CuPy-based simple linear regression for many response vectors (rows of `y`)
	against a single predictor `x`, allowing for NaNs in `y`.  Any sample whose
	response is NaN is ignored for that region; its x-value is excluded from all
	sums for that region as well.

	Parameters
	----------
	x : (n_samples,) cupy.ndarray
	    Predictor values.
	y : (n_regions, n_samples) cupy.ndarray
	    Each row is the response for one region.

	Returns
	-------
	pvalue, slope : (n_regions,) cupy.ndarray
	'''
	mask = ~cp.isnan(y)
	n_eff = mask.sum(axis=1)

	valid = n_eff >= 3
	n_eff = n_eff.astype(cp.float32)

	xb = x[None, :]

	sum_x = (xb * mask).sum(axis=1)
	sum_y = cp.where(mask, y, 0).sum(axis=1)
	mean_x = cp.where(valid, sum_x / n_eff, cp.nan)
	mean_y = cp.where(valid, sum_y / n_eff, cp.nan)

	x_dev = (xb - mean_x[:, None]) * mask
	y_dev = (y  - mean_y[:, None]) * mask

	ss_xx = (x_dev ** 2).sum(axis=1)
	ss_xy = (x_dev * y_dev).sum(axis=1)
	ss_yy = (y_dev ** 2).sum(axis=1)

	slope = ss_xy / ss_xx

	r = ss_xy / cp.sqrt(ss_xx * ss_yy)
	r = cp.clip(r, -1.0, 1.0)

	t_stat = r * cp.sqrt((n_eff - 2) / (1.0 - r ** 2))
	v = n_eff - 2.0
	pvalue = cpx_spl.betainc(v / 2.0, 0.5, v / (v + t_stat ** 2))

	missing = (~valid) | cp.isnan(r) | cp.isnan(pvalue)
	slope[missing] = 0
	pvalue[missing] = 1
	return pvalue, slope


def compute_chromosome_gene_regressions(het_df, total_snp_df, sample_dates, n_iterations = 1001, random_state = None):
	'''
	Run CuPy-accelerated linear regressions of per-sample heterozygosity
	on sample collection dates for every genomic element, with bootstrap
	resampling.
	'''

	het_counts   = het_df.iloc[:, 4:].to_numpy(dtype = np.float32)
	total_counts = total_snp_df.iloc[:, 4:].to_numpy(dtype = np.float32)
	het_matrix   = het_counts / total_counts

	n_regions, n_samples = het_matrix.shape

	if random_state is not None:
		cp.random.seed(random_state)

	bootstrap = cp.array(cp.random.choice(
		n_samples, (n_iterations, n_samples), replace = True
	), dtype=cp.int32)
	bootstrap[0] = cp.arange(n_samples, dtype = cp.int32)

	regressions = np.empty((n_iterations, n_regions, 2), dtype = np.float32)

	het_gpu   = cp.asarray(het_matrix)
	dates_gpu = cp.asarray(sample_dates)

	for i in range(n_iterations):
		idx = bootstrap[i]
		x   = dates_gpu[idx]
		y   = het_gpu[:, idx]

		pvals, slopes = cupy_parallel_linregress(x, y)

		regressions[i, :, 0] = cp.asnumpy(slopes)
		regressions[i, :, 1] = cp.asnumpy(pvals)

	return regressions


def benjamini_hochberg_mask(pvals, alpha = 0.01):
	'''
	Return boolean mask of the positions that survive Benjamini-Hochberg correction.
	'''
	pvals = np.asarray(pvals)
	shape = pvals.shape
	if pvals.size == 0:
		return np.zeros(shape, dtype=bool)

	flat = pvals.ravel()
	m = flat.size
	sort_idx = np.argsort(flat)
	sorted_p = flat[sort_idx]
	thresholds = alpha * (np.arange(1, m + 1) / m)
	ok = sorted_p <= thresholds

	if not ok.any():
		return np.zeros(shape, dtype=bool)

	k_max = np.where(ok)[0].max()
	reject_flat = np.zeros(m, dtype=bool)
	reject_flat[sort_idx[:k_max + 1]] = True
	return reject_flat.reshape(shape)


def save_regression_results(regressions, regions, output_path, fdr = 0.01):
	regressions = cp.asnumpy(regressions)
	slopes = regressions[..., 0]
	pvals  = regressions[..., 1]
	sig_primary   = benjamini_hochberg_mask(pvals[0], alpha=fdr)
	sig_bootstrap = benjamini_hochberg_mask(pvals[1:], alpha=fdr).sum(axis=0)
	result_df = pd.concat(
		[regions.reset_index(drop = True),
		 pd.DataFrame({
			'slope'            : slopes[0],
			'beta_H'           : -1e3 * slopes[0],
			'p_value'          : pvals[0],
			'significant'      : sig_primary,
			'bootstrap_repeats': sig_bootstrap
		 })],
		axis = 1
	)
	os.makedirs(output_path, exist_ok=True)
	result_df.to_csv(f'{output_path}/regressions.csv', sep='\t', index=False)


def process_genotypes_to_heterozygosity(args):
	labelled_regions = parse_labelled_regions(args.labelled_regions)

	samples = pd.read_csv(args.samples_csv, sep = '\t', index_col = 0)
	sample_dates = samples['date'].to_numpy()

	regions = pd.read_csv(
		args.regions_file, sep = '\t',
		names = ['chr', 'start', 'end', 'name'],
		dtype = {'chr': str, 'start': int, 'end': int, 'name': str}
	)

	bed_file, sample_bed_indices, bim = init_bed_reader(
		args.bed_prefix, samples.index, args.threads
	)

	het_win_df, total_win_df = compute_window_het(
		bed_file, sample_bed_indices, bim, samples.index,
		window_size = args.window_size,
		step_size = args.window_size,
		exclude_regions = [labelled_regions['mhc']]
	)

	het_win_df.to_csv(f'{args.output_path}/windows_het_counts.csv',  sep='\t', index=False)
	total_win_df.to_csv(f'{args.output_path}/windows_total_counts.csv', sep='\t', index=False)

	wg_het = het_win_df.iloc[:, 3:].sum(0) / total_win_df.iloc[:, 3:].sum(0)
	wg_het.name = 'wg_het'

	het_lab_cnt_df, tot_lab_cnt_df = compute_labelled_counts(
		bed_file, sample_bed_indices, bim, samples.index, labelled_regions
	)
	labelled_frac_df = compute_labelled_het(het_lab_cnt_df, tot_lab_cnt_df, samples.index)

	samples_het = pd.concat([samples, labelled_frac_df, wg_het], axis = 1)
	os.makedirs(args.output_path, exist_ok = True)
	samples_het.to_csv(f'{args.output_path}/samples_het.csv', sep = '\t')

	het_df, total_snp_df = compute_het(
		bed_file, sample_bed_indices, bim, samples.index, regions
	)
	regressions = compute_chromosome_gene_regressions(
		het_df, total_snp_df, sample_dates
	)
	save_regression_results(regressions, regions, args.output_path)

	het_df.to_csv(f'{args.output_path}/regions_het_snps.csv', sep = '\t', index = False)
	total_snp_df.to_csv(f'{args.output_path}/regions_total_snps.csv', sep = '\t', index = False)

	null_long_df = selection_coeff_null(
		het_win_df, total_win_df, sample_dates,
		time_window_start = args.time_window_start,
		time_window    = args.time_window,
		time_step      = args.time_step
	)
	null_long_df.to_csv(f'{args.output_path}/null_selection_coeff.csv', sep = '\t', index = False)

	label_long_df = selection_coeff_null(
		het_lab_cnt_df, tot_lab_cnt_df, sample_dates,
		time_window_start = args.time_window_start,
		time_window    = args.time_window,
		time_step      = args.time_step
	)
	label_long_df.to_csv(f'{args.output_path}/labelled_selection_coeff.csv', sep = '\t', index = False)


def main():
	parser = argparse.ArgumentParser(description="Compute heterozygosity for genes/windows")
	parser.add_argument('--samples-csv', type=str, required=True, help='Sample metadata corresponding to data in the PLINK files (all samples should be in the FAM file with the same IDs)')
	parser.add_argument('--bed-prefix', type=str, required=True, help='Prefix of PLINK BED/FAM/BIM files containing genotypes')
	parser.add_argument('--labelled-regions', type=str, required=True, help='List of large genomic regions (as opposed to fixed windows or genes) to add as columns to samples table and compute selection coefficients for')
	parser.add_argument('--regions-file', type=str, required=True, help='TSV with list of genomic regions, columns are <chr> <start> <end> <name> (BED format)')
	parser.add_argument('--window-size', type = int, default=1000000, help='Genomic window size for computing null distribution of selection coefficients')
	parser.add_argument('--time-window-start', type = int, default=12000, help='First temporal window for selection coefficient analysis')
	parser.add_argument('--time-window', type = int, default=2000, help='Size of temporal window for computing mean heterozygosity for selection coefficient analysis')
	parser.add_argument('--time-step', type = int, default=1000, help='Distance between two adjacent windows of size `time_window`')
	parser.add_argument('--threads', type=int, default=4)
	parser.add_argument('--output-path', type=str, required=True)
	args = parser.parse_args()
	process_genotypes_to_heterozygosity(args)


if __name__ == '__main__':
	main()
