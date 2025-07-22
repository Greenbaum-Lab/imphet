import subprocess
from pathlib import Path


def _run(cmd):
	'''Run a command and raise if it fails.'''
	completed = subprocess.run(cmd, shell=isinstance(cmd, str), text=True)
	if completed.returncode:
		raise RuntimeError(f'command failed ({completed.returncode}): {cmd}')


def _index_bcf(bcf_path, threads=1):
	if not bcf_path.with_suffix(bcf_path.suffix + '.csi').exists():
		_run(['bcftools', 'index', '-f', '--threads', str(threads), str(bcf_path)])


def concat_sample(sample_dir):
	'''
	Concatenate chr1-22 per-sample BCFs into one file and index it.
	Returns the path to the concatenated BCF.
	'''
	sample_id = sample_dir.name
	output_bcf = sample_dir / f'{sample_id}.bcf'
	if output_bcf.exists():
		return output_bcf

	input_paths = [
		str(sample_dir / f'{sample_id}_{c}.bcf')
		for c in range(1, 23)
		if (sample_dir / f'{sample_id}_{c}.bcf').exists()
	]

	_run(['bcftools', 'concat', '--threads', '1', '-Ob', '-o', str(output_bcf), *input_paths])
	_index_bcf(output_bcf)
	return output_bcf


def merge_samples(sample_dirs, prefix, threads=32):
	'''
	Merge all per-sample BCFs into one multi-sample BCF.
	Returns the merged BCF path.
	'''
	output_bcf = prefix.with_suffix('.bcf')
	if output_bcf.exists():
		return output_bcf

	input_paths = [str(directory / f'{directory.name}.bcf') for directory in sample_dirs]
	_run([
		'bcftools', 'merge',
		'--threads', str(threads),
		'--force-samples',
		'-Ob', '-o', str(output_bcf),
		*input_paths
	])
	_index_bcf(output_bcf, threads)
	return output_bcf


def setgt_filter(in_bcf, threads=32, min_gp=None, min_dp=None, expr=None):
	'''
	Blank genotypes with bcftools +setGT.
	If expr is provided it is used as-is.
	Otherwise the first non-None of min_gp or min_dp is used to build the expression:
	    min_gp -> "SMPL_MAX(FMT/GP)<min_gp"
	    min_dp -> "FMT/DP<min_dp"
	Returns the filtered BCF path. If no filtering is requested the input is returned.
	'''
	if expr is None:
		if min_gp is not None:
			expr = f'SMPL_MAX(FMT/GP)<{min_gp}'
		elif min_dp is not None:
			expr = f'FMT/DP<{min_dp}'

	if expr is None:
		return in_bcf

	filtered_bcf = in_bcf.with_name(in_bcf.stem + '_filtered.bcf')
	if filtered_bcf.exists():
		return filtered_bcf

	cmd = f'bcftools +setGT {in_bcf} -Ob -o {filtered_bcf} --threads {threads} -- -t q -n . -i "{expr}"'
	_run(cmd)
	_index_bcf(filtered_bcf, threads)
	return filtered_bcf


def plink_convert(in_bcf, prefix, threads=32):
	'''
	Convert a BCF to PLINK binary format. Skips if *.bed exists.
	'''
	bed_path = prefix.with_suffix('.bed')
	if bed_path.exists():
		return

	mem_mb = threads * 7000
	_run([
		'plink',
		'--bcf', str(in_bcf),
		'--threads', str(threads),
		'--memory', str(mem_mb),
		'--make-bed',
		'--allow-no-sex',
		'--keep-allele-order',
		'--indiv-sort', '0',
		'--out', str(prefix)
	])
