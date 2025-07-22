#!/usr/bin/env python3

import argparse
import subprocess
import numpy as np
import pandas as pd
import pysam
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def run_cmd(cmd):
	'''Execute a shell or list command and raise on failure.'''
	subprocess.run(cmd, shell=isinstance(cmd, str), check=True, text=True)


def load_coverage(cov_path):
	'''Return coverage vector from a BED file produced by bedtools map.'''
	return (
		pd.read_csv(
			cov_path,
			sep='\t',
			header=None,
			usecols=[3],
			names=['cov'],
			na_values=['.']
		)
		.fillna(0)
		['cov']
		.values
	)


def generate_coverage(bam_path, ref_fai, cache_dir, chrom_label, bin_size):
	'''
	Create binned coverage BED for bam_path if absent and return its Path.
	'''
	sample_id = bam_path.stem
	out_path = cache_dir / f'{sample_id}.{chrom_label}.{bin_size}.cov.bed'
	if out_path.exists():
		return out_path

	tmp_per_base = out_path.with_suffix('.per_base.bed')
	bins_file = cache_dir / f'bins_{chrom_label}_{bin_size}.bed'

	if not bins_file.exists():
		run_cmd(
			f'bedtools makewindows -g {ref_fai} -w {bin_size} | '
			f"awk '$1==\"{chrom_label}\"' > {bins_file}"
		)

	run_cmd(
		f'samtools view -b {bam_path} {chrom_label} | '
		f'bedtools genomecov -ibam stdin -bg | '
		f'sort -k1,1 -k2,2n > {tmp_per_base}'
	)

	run_cmd(['bedtools', 'map', '-a', str(bins_file), '-b', str(tmp_per_base),
			 '-c', '4', '-o', 'mean', '>', str(out_path)])
	tmp_per_base.unlink(missing_ok=True)
	return out_path


def downsample_bam(origin_bam, template_bam, coverage_dir, chrom_label, bin_size, out_bam):
	'''
	Write a down-sampled BAM derived from origin_bam to out_bam so that its
	coverage follows template_bam on chrom_label.
	'''
	origin_cov = load_coverage(
		coverage_dir / f'{origin_bam.stem}.{chrom_label}.{bin_size}.cov.bed'
	)
	template_cov = load_coverage(
		coverage_dir / f'{template_bam.stem}.{chrom_label}.{bin_size}.cov.bed'
	)
	proportions = np.where(origin_cov > 0, template_cov / origin_cov, 0.0)

	bam_in = pysam.AlignmentFile(origin_bam, 'rb')
	bam_out = pysam.AlignmentFile(out_bam, 'wb', template=bam_in)

	for read in bam_in.fetch(chrom_label):
		start = read.reference_start
		if start < 0:
			continue
		bin_idx = start // bin_size
		if bin_idx >= len(proportions):
			continue
		if np.random.rand() <= proportions[bin_idx]:
			bam_out.write(read)

	bam_in.close()
	bam_out.close()
	run_cmd(['samtools', 'index', '-b', str(out_bam)])


def build_coverage_cache(bam_paths, ref_fai, cache_dir, chrom_label, bin_size, threads):
	'''
	Ensure coverage files exist for every BAM in bam_paths.
	'''
	with ThreadPoolExecutor(max_workers=threads) as pool:
		futures = [
			pool.submit(
				generate_coverage,
				bam_file,
				ref_fai,
				cache_dir,
				chrom_label,
				bin_size
			)
			for bam_file in bam_paths
		]
		for fut in as_completed(futures):
			fut.result()


def downsample_batch(origin_dir, template_dir, coverage_dir, ref_fai,
					 chrom_label, out_dir, bin_size=50, threads=4):
	'''
	Generate down-sampled BAMs for one chromosome and return their paths.
	'''
	origin_bams = sorted(Path(origin_dir).glob(f'*{chrom_label}.bam'))
	template_bams = sorted(Path(template_dir).glob(f'*{chrom_label}.bam'))
	if not origin_bams or not template_bams:
		raise FileNotFoundError('No BAMs found for requested chromosome')

	coverage_dir = Path(coverage_dir)
	coverage_dir.mkdir(parents=True, exist_ok=True)

	build_coverage_cache(origin_bams + template_bams, ref_fai, coverage_dir,
						 chrom_label, bin_size, threads)

	out_dir = Path(out_dir)
	out_dir.mkdir(parents=True, exist_ok=True)

	out_paths = []
	with ThreadPoolExecutor(max_workers=threads) as pool:
		futures = []
		for orig_bam in origin_bams:
			for tmpl_bam in template_bams:
				out_name = f'{orig_bam.stem}_{tmpl_bam.stem}_{chrom_label}.bam'
				out_path = out_dir / out_name
				if out_path.exists():
					out_paths.append(out_path)
					continue
				out_paths.append(out_path)
				futures.append(
					pool.submit(
						downsample_bam,
						orig_bam,
						tmpl_bam,
						coverage_dir,
						chrom_label,
						bin_size,
						out_path
					)
				)
		for fut in as_completed(futures):
			fut.result()
	return out_paths


def parse_args():
	parser = argparse.ArgumentParser(description='Downsample BAMs by coverage.')
	parser.add_argument('--origin-dir', required=True)
	parser.add_argument('--template-dir', required=True)
	parser.add_argument('--coverage-dir', required=True)
	parser.add_argument('--reference-fai', required=True)
	parser.add_argument('--chromosome', required=True)
	parser.add_argument('--output-dir', required=True)
	parser.add_argument('--bin-size', type=int, default=50)
	parser.add_argument('--threads', type=int, default=4)
	return parser.parse_args()


def main():
	args = parse_args()
	downsample_batch(
		args.origin_dir,
		args.template_dir,
		args.coverage_dir,
		args.reference_fai,
		args.chromosome,
		args.output_dir,
		args.bin_size,
		args.threads
	)


if __name__ == '__main__':
	main()
