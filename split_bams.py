import argparse
import os
import sys
import subprocess
import tempfile
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed


def index_bam(bam_file, threads):
	'''
	Create an index for a BAM if it does not already exist.
	'''
	bai_file = f'{bam_file}.bai'
	if not os.path.exists(bai_file):
		print(f'[INFO] Indexing {bam_file}')
		sys.stdout.flush()
		subprocess.run(['samtools', 'index', '-@', str(threads), bam_file], check=True)


def sort_bam(input_bam, output_bam, threads):
	'''
	Sort a BAM file.
	'''
	print(f'[INFO] Sorting {input_bam} -> {output_bam}')
	sys.stdout.flush()
	subprocess.run(
		['samtools', 'sort', '-@', str(threads), '-o', output_bam, input_bam],
		check=True
	)


def extract_chromosome(source_bam, chromosome, output_dir, genetic_id, threads):
	'''
	Extract reads for one chromosome and create an index for the result.
	'''
	chrom_dir = os.path.join(output_dir, chromosome)
	os.makedirs(chrom_dir, exist_ok=True)

	out_bam = os.path.join(chrom_dir, f'{genetic_id}.bam')

	print(f'[INFO] Extracting {chromosome} for {genetic_id} -> {out_bam}')
	sys.stdout.flush()

	subprocess.run(
		[
			'samtools', 'view',
			'-@', str(threads),
			'-b',
			'-o', out_bam,
			source_bam,
			chromosome.replace('chr', '')
		],
		check=True
	)

	index_bam(out_bam, threads)


def group_samples(input_dir, samples_file):
	'''
	Build a dict mapping genetic_id to list of BAM paths.
	'''
	samples = pd.read_csv(samples_file, delimiter='\t', header=0)
	samples['bam_path'] = samples['accession'].apply(
		lambda acc: os.path.join(input_dir, f'{acc}.bam')
	)
	samples = samples[samples['bam_path'].apply(os.path.exists)]
	return (
		samples.groupby('genetic_id')['bam_path']
		.apply(list)
		.to_dict()
	)


def process_sample(sample_id, bam_files, chromosomes, output_dir, threads):
	is_merged = len(bam_files) > 1

	with tempfile.TemporaryDirectory() as temp_dir:
		if is_merged:
			merged_bam = os.path.join(temp_dir, f'{sample_id}.merged.bam')
			print(f'[INFO] Merging accessions for {sample_id} -> {merged_bam}')
			sys.stdout.flush()
			subprocess.run(
				['samtools', 'merge', '-@', str(threads), merged_bam] + bam_files,
				check=True
			)
			sorted_bam = os.path.join(temp_dir, f'{sample_id}.sorted.bam')
			sort_bam(merged_bam, sorted_bam, threads * len(chromosomes))
			index_bam(sorted_bam, threads * len(chromosomes))
			source_bam = sorted_bam
		else:
			source_bam = bam_files[0]
			index_bam(source_bam, threads * len(chromosomes))
			print(f'[INFO] Single accession for {sample_id}: {source_bam}')
			sys.stdout.flush()

		with ThreadPoolExecutor(max_workers=len(chromosomes)) as executor:
			futures = [
				executor.submit(
					extract_chromosome,
					source_bam,
					chrom,
					output_dir,
					sample_id,
					threads
				)
				for chrom in chromosomes
			]
			for f in as_completed(futures):
				f.result()

		for chrom in chromosomes:
			bam_path = os.path.join(output_dir, chrom, f'{sample_id}.bam')
			bai_path = f'{bam_path}.bai'
			if not (os.path.exists(bam_path) and os.path.exists(bai_path)):
				raise RuntimeError(f'Missing output for {sample_id} {chrom}')

		for bam in bam_files:
			if os.path.exists(bam):
				os.remove(bam)
			bai = f'{bam}.bai'
			if os.path.exists(bai):
				os.remove(bai)


def main():
	parser = argparse.ArgumentParser(
		description='Merge accessions, split by chromosome, and clean up.'
	)
	parser.add_argument('-i', '--input-dir', required=True, help='Directory containing BAM files')
	parser.add_argument('-o', '--output-dir', default='.', help='Directory to store output BAM files')
	parser.add_argument('-s', '--samples-file', required=True, help='Tab-delimited sample file with genetic IDs')
	parser.add_argument('-t', '--threads', type=int, default=1, help='Threads per samtools call (per chromosomes); assumes you have `len(chromosomes) * threads` overall')
	parser.add_argument(
		'-c',
		'--chromosomes',
		nargs='+',
		default=[f'{x}' for x in range(1, 23)],
		help='Chromosomes to process'
	)

	args = parser.parse_args()

	os.makedirs(args.output_dir, exist_ok=True)

	sample_groups = group_samples(args.input_dir, args.samples_file)

	for sample_id, bam_files in sample_groups.items():
		print(f'[INFO] Processing sample {sample_id}')
		sys.stdout.flush()
		try:
			process_sample(
				sample_id,
				bam_files,
				args.chromosomes,
				args.output_dir,
				args.threads
			)
		except Exception as e:
			print(f'[ERROR] Sample {sample_id} failed: {e}')
			sys.stdout.flush()
			continue

	print('[INFO] All samples processed')
	sys.stdout.flush()


if __name__ == '__main__':
	main()
