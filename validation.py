#!/usr/bin/env python3
'''
Validation pipeline

High-coverage source sample workflow
	1. Down-sample BAMs to template coverage
	2. Call genotypes on down-sampled BAMs and depth-filter
	3. Call genotypes on original BAMs for truth set and depth-filter
	4. Impute the down-sampled BAMs with GLIMPSE2 and probability-filter
	5. Merge imputed and truth BCFs, convert to PLINK
'''

import argparse
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from downsample import downsample_batch
from glimpse_wrapper import prepare_reference
from impute_genotypes import impute_genotypes
from bcfs_to_plink import (
	concat_sample,
	plink_convert,
	setgt_filter,
	merge_samples,
)
from utils import create_bcf_from_bam


def _run(cmd):
	subprocess.run(cmd, shell=isinstance(cmd, str), check=True)


def prepare_references(chrom_labels, glimpse_root, reference_prefix):
	for chrom in chrom_labels:
		prepare_reference(glimpse_root, reference_prefix, chrom)


def downsample_all(source_dir, template_dir, ref_fai, chrom_labels,
				   dst_root, bin_size, threads):
	down_root = dst_root / 'downsampled'
	down_root.mkdir(parents=True, exist_ok=True)
	down_bam_paths = []
	for chrom in chrom_labels:
		chrom_dst = down_root / f'chr{chrom}'
		chrom_dst.mkdir(exist_ok=True)
		down_bam_paths += downsample_batch(
			origin_dir=Path(source_dir),
			template_dir=Path(template_dir),
			coverage_dir=dst_root / 'cov',
			ref_fai=ref_fai,
			chrom_label=f'chr{chrom}',
			out_dir=chrom_dst,
			bin_size=bin_size,
			threads=threads
		)
	return down_root, down_bam_paths


def _sites_file(reference_prefix, chromosome):
	prefix_name = Path(reference_prefix).name
	return (
		Path(reference_prefix).parent /
		f'{prefix_name}_{chromosome}' /
		f'{prefix_name}_{chromosome}.sites.vcf.gz'
	)


def call_genotypes(bam_paths, reference_prefix, out_dir, threads):
	out_dir.mkdir(parents=True, exist_ok=True)

	def worker(bam_path):
		chromosome = bam_path.parent.name[3:]
		ref_sites = _sites_file(reference_prefix, chromosome)
		sample_id = bam_path.stem.split('_')[0]
		sample_dir = out_dir / sample_id
		sample_dir.mkdir(exist_ok=True)
		out_bcf = sample_dir / f'{sample_id}_{chromosome}.bcf'
		if out_bcf.exists():
			return
		create_bcf_from_bam(
			input_bam=str(bam_path),
			chromosome=chromosome,
			reference_fasta=str(reference_prefix) + '.fasta',  # .fasta path derived from prefix
			reference_sites=str(ref_sites),
			output_path=str(out_bcf),
			tsv=None
		)

	with ThreadPoolExecutor(max_workers=threads) as pool:
		list(pool.map(worker, bam_paths))


def prepare_called_route(called_root, depth_thr, threads):
	for sample_dir in called_root.iterdir():
		if not sample_dir.is_dir():
			continue
		concat_sample(sample_dir)
		bcf_path = sample_dir / f'{sample_dir.name}.bcf'
		filtered = setgt_filter(bcf_path, threads=threads, min_dp=depth_thr)
		plink_dir = called_root.parent / 'plink' / 'called'
		plink_dir.mkdir(parents=True, exist_ok=True)
		plink_convert(filtered, plink_dir / sample_dir.name, threads)


def create_truth_bcfs(source_bams, reference_prefix, truth_root,
					  depth_thr, threads):
	call_genotypes(source_bams, reference_prefix, truth_root, threads)

	for sample_dir in truth_root.iterdir():
		if not sample_dir.is_dir():
			continue
		concat_sample(sample_dir)
		bcf_path = sample_dir / f'{sample_dir.name}.bcf'
		filtered = setgt_filter(bcf_path, threads=threads, min_dp=depth_thr)
		shutil.move(filtered, bcf_path.with_name(filtered.name))
		shutil.move(filtered.with_suffix('.bcf.csi'),
					bcf_path.with_name(filtered.stem + '.csi'))


def process_imputed(imputed_root, prob_thr, threads):
	for sample_dir in imputed_root.iterdir():
		if not sample_dir.is_dir():
			continue
		concat_bcf = sample_dir / f'{sample_dir.name}.bcf'
		if concat_bcf.exists():
			filtered = setgt_filter(concat_bcf, threads=threads, min_gp=prob_thr)
			plink_dir = imputed_root.parent / 'plink' / 'imputed_separate'
			plink_dir.mkdir(parents=True, exist_ok=True)
			plink_convert(filtered, plink_dir / sample_dir.name, threads)


def merge_imputed_with_truth(imputed_root, truth_root, output_root,
							 prob_thr, threads):
	merge_dirs = [d for d in imputed_root.iterdir() if d.is_dir()]
	merge_dirs += [d for d in truth_root.iterdir() if d.is_dir()]
	merged_prefix = imputed_root / 'merged_full'
	merged_bcf = merge_samples(merge_dirs, merged_prefix, threads)
	filtered = setgt_filter(merged_bcf, threads=threads, min_gp=prob_thr)
	final_plink_dir = output_root / 'plink' / 'imputed'
	final_plink_dir.mkdir(parents=True, exist_ok=True)
	plink_convert(filtered, final_plink_dir / 'merged_full', threads)


def run_validation_pipeline(src_dir, tmpl_dir, reference_prefix, glimpse_root,
							ref_fasta, out_dir, chromosome, bin_sz,
							prob_thr, depth_thr, threads):
	tmp_root = Path(tempfile.gettempdir()) / 'glimpse_val'
	tmp_root.mkdir(exist_ok=True)

	src_bams = [p for p in Path(src_dir).rglob('*.bam') if p.parent.name.startswith('chr')]
	chrom_set = {p.parent.name[3:] for p in src_bams}

	prepare_references(chrom_set, glimpse_root, reference_prefix)

	down_root, down_bams = downsample_all(src_dir, tmpl_dir,
										  Path(ref_fasta).with_suffix('.fai'),
										  chrom_set, out_dir,
										  bin_sz, threads)

	called_root = out_dir / 'called'
	call_genotypes(down_bams, reference_prefix, called_root, threads)

	prepare_called_route(called_root, depth_thr, threads)

	imputed_root = out_dir / 'imputed'
	impute_genotypes(down_root, reference_prefix, imputed_root,
					 glimpse_root, threads, prob_thr)

	process_imputed(imputed_root, prob_thr, threads)

	truth_root = out_dir / 'truth'
	truth_root.mkdir(parents=True, exist_ok=True)
	create_truth_bcfs(src_bams, reference_prefix, truth_root,
					  depth_thr, threads)

	merge_imputed_with_truth(imputed_root, truth_root, out_dir,
							 prob_thr, threads)

	print(f'Validation finished, results in {out_dir}')


def main():
	parser = argparse.ArgumentParser(description='Validation pipeline for GLIMPSE2')
	parser.add_argument('--source-bam-directory', required=True)
	parser.add_argument('--template-bam-directory', required=True)
	parser.add_argument('--reference-prefix', required=True)
	parser.add_argument('--glimpse-directory', required=True)
	parser.add_argument('--reference-fasta', required=True)
	parser.add_argument('--output-dir', required=True)
	parser.add_argument('--chromosome', type=int, required=True)
	parser.add_argument('--bin-size', type=int, default=50)
	parser.add_argument('--threads', type=int, default=8)
	parser.add_argument('--min-gp', type=float, default=0.90)
	parser.add_argument('--min-dp', type=int, default=4)
	args = parser.parse_args()

	run_validation_pipeline(
		Path(args.source_bam_directory).resolve(),
		Path(args.template_bam_directory).resolve(),
		Path(args.reference_prefix).resolve(),
		Path(args.glimpse_directory).resolve(),
		Path(args.reference_fasta).resolve(),
		Path(args.output_dir).resolve(),
		args.chromosome,
		args.bin_size,
		args.min_gp,
		args.min_dp,
		args.threads
	)


if __name__ == '__main__':
	main()
