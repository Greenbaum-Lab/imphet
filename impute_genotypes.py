#!/usr/bin/env python3
'''
End-to-end GLIMPSE2 imputation pipeline.

1. Prepares each chromosome of the reference once,
2. Imputes every BAM found in <source_bam_dir>/chr*/,
3. Concatenates the 22 chromosome-BCFs per sample,
4. Merges all samples, applies a +setGT filter,
5. Converts the final BCF to PLINK.

Output layout (under --output-dir):

	<output_dir>/<sample_id>/<sample_id>_<chr>.bcf
	<output_dir>/<sample_id>/<sample_id>.bcf
	<output_dir>/merged_full.bcf
	<output_dir>/merged_full_filtered.bcf
	<output_dir>/merged_full.bed/.bim/.fam
'''

import argparse
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from glimpse_wrapper import prepare_reference, impute_sample
from bcfs_to_plink import (
	concat_sample,
	merge_samples,
	setgt_filter,
	plink_convert,
)


def _discover_bams(bam_root):
	bams = []
	for chr_dir in bam_root.iterdir():
		if chr_dir.is_dir() and chr_dir.name.startswith('chr'):
			bams.extend(chr_dir.glob('*.bam'))
	return bams


def _impute_all(bam_root, ref_prefix, out_root, glimpse_dir, threads):
	out_root.mkdir(parents=True, exist_ok=True)
	tmp_root = Path(tempfile.gettempdir()) / 'glimpse_work'
	tmp_root.mkdir(parents=True, exist_ok=True)

	bam_paths = _discover_bams(bam_root)

	chromosomes = {p.parent.name[3:] for p in bam_paths}
	for chrom in sorted(chromosomes):
		prepare_reference(glimpse_dir, ref_prefix, chrom)

	def _worker(bam_path):
		sample_id = bam_path.stem
		chrom = bam_path.parent.name[3:]
		target_bcf = out_root / sample_id / f'{sample_id}_{chrom}.bcf'
		if target_bcf.exists() and target_bcf.with_suffix(target_bcf.suffix + '.csi').exists():
			print(f'skip {sample_id} chr{chrom} – already imputed')
			return
		try:
			final_bcf = impute_sample(glimpse_dir, ref_prefix, chrom, bam_path, tmp_root)

			sample_dir = out_root / sample_id
			sample_dir.mkdir(exist_ok=True)

			shutil.move(final_bcf, target_bcf)

			idx_src = final_bcf.with_suffix(final_bcf.suffix + '.csi')
			if idx_src.exists():
				shutil.move(idx_src, target_bcf.with_suffix(target_bcf.suffix + '.csi'))

			#shutil.rmtree(final_bcf.parent, ignore_errors=True)
			print(f'finished {sample_id} chr{chrom}')
		except Exception as exc:
			print(f'sample {sample_id} chr{chrom} failed: {exc}')
			raise

	with ThreadPoolExecutor(max_workers=threads) as pool:
		list(pool.map(_worker, bam_paths))


def _postprocess(out_root, threads, min_gp):
	sample_dirs = [d for d in out_root.iterdir() if d.is_dir()]

	for sample_dir in sample_dirs:
		concat_sample(sample_dir)

	merged_prefix = out_root / 'merged_full'
	merged_bcf = merge_samples(sample_dirs, merged_prefix, threads)
	filtered_bcf = setgt_filter(merged_bcf, threads=threads, min_gp=min_gp)
	plink_convert(filtered_bcf, merged_prefix, threads=threads)


def impute_genotypes(
	source_bam_dir,
	reference_dir,
	output_dir,
	glimpse_dir,
	threads,
	min_gp,
):
	_impute_all(
		source_bam_dir,
		reference_dir,
		output_dir,
		glimpse_dir,
		threads,
	)
	_postprocess(output_dir, threads, min_gp)


def main():
	p = argparse.ArgumentParser(description='Impute genotypes with GLIMPSE2')
	p.add_argument('--source-bam-directory', required=True)
	p.add_argument('--reference-prefix', required=True)
	p.add_argument('--output-dir', required=True)
	p.add_argument('--glimpse-directory', required=True)
	p.add_argument('--threads', type=int, default=8)
	p.add_argument('--min-gp', type=float, default=0.90)
	args = p.parse_args()

	impute_genotypes(
		Path(args.source_bam_directory).resolve(),
		Path(args.reference_prefix).resolve(),
		Path(args.output_dir).resolve(),
		Path(args.glimpse_directory).resolve(),
		args.threads,
		args.min_gp,
	)


if __name__ == '__main__':
	main()
