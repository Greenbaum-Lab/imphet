import argparse
import os
import sys
import subprocess
import tempfile
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

def index_bam(bam_file):
	'''
	Index the BAM file if the corresponding BAI file doesn't exist.
	'''
	bai_file = f'{bam_file}.bai'
	if not os.path.exists(bai_file):
		print(f'[INFO] Indexing: {bam_file}')
		sys.stdout.flush()
		cmd = ['samtools', 'index', bam_file]
		subprocess.run(cmd, check=True)
	else:
		print(f'[INFO] BAI found, skipping indexing: {bam_file}')
		sys.stdout.flush()

def extract_chromosome(bam_file, chromosome, output_dir, genetic_id, threads = 1):
	'''
	Extracts reads for a given chromosome from a BAM file into a separate BAM.
	'''
	chrom_output_dir = os.path.join(output_dir, chromosome)
	os.makedirs(chrom_output_dir, exist_ok=True)

	out_bam = os.path.join(chrom_output_dir, f'{genetic_id}.bam')
	cmd = [
		'samtools', 'view',
		'-@', str(threads),
		'-b',
		'-o', out_bam,
		bam_file,
		chromosome.replace('chr', '')
	]

	print(f'[INFO] Extracting {chromosome} for {genetic_id} to {out_bam}')
	sys.stdout.flush()
	subprocess.run(cmd, check=True)

def sort_and_index_bam_if_needed(bam_file, threads, is_merged, temp_dir):
	'''
	Sorts and indexes a BAM file only if it was merged.
	'''
	if not is_merged:
		print(f'[INFO] BAM file is already sorted and indexed: {bam_file}')
		sys.stdout.flush()
		return bam_file

	sorted_bam = os.path.join(temp_dir, os.path.basename(bam_file).replace('.bam', '.sorted.bam'))
	print(f'[INFO] Sorting BAM file: {bam_file} -> {sorted_bam}')
	sys.stdout.flush()
	subprocess.run(['samtools', 'sort', '-@', str(threads), '-o', sorted_bam, bam_file], check=True)

	print(f'[INFO] Indexing sorted BAM file: {sorted_bam}')
	sys.stdout.flush()
	subprocess.run(['samtools', 'index', sorted_bam], check=True)

	return sorted_bam

def get_bam_files(input_dir):
	'''
	Retrieve a list of BAM files from the specified directory.
	'''
	return [
		os.path.join(input_dir, file)
		for file in os.listdir(input_dir)
		if file.endswith('.bam')
	]

def group_samples(input_dir, samples_file, bam_files):
	'''
	Load sample information and group BAM files by genetic ID.
	'''
	samples = pd.read_csv(samples_file, delimiter='\t', header=0)
	grouped = {}
	for i, sample in samples.iterrows():
		bam_filename = f'{input_dir}/{sample.accession}.bam'
		file_exists = bam_filename in bam_files
		if file_exists:
			if not sample.genetic_id in grouped:
				grouped[sample.genetic_id] = []
			grouped[sample.genetic_id].append(bam_filename)
	return grouped

def process_sample(sample_id, bam_files, chromosomes, output_dir, threads):
	is_merged = len(bam_files) > 1

	with tempfile.TemporaryDirectory() as temp_dir:
		if not is_merged:
			combined_bam = bam_files[0]
			print(f'[INFO] Single accession for sample {sample_id}, skipping merge: {combined_bam}')
			sys.stdout.flush()
		else:
			combined_bam = os.path.join(temp_dir, f'{sample_id}.bam')
			cmd_merge = ['samtools', 'merge', '-@', str(threads), combined_bam] + bam_files
			print(f'[INFO] Merging accessions for sample {sample_id} into {combined_bam}')
			sys.stdout.flush()
			subprocess.run(cmd_merge, check=True)

		sorted_bam = sort_and_index_bam_if_needed(combined_bam, threads, is_merged, temp_dir)

		with ThreadPoolExecutor(max_workers=len(chromosomes)) as executor:
			future_to_chromosome = {
				executor.submit(
					extract_chromosome,
					sorted_bam,
					chrom,
					output_dir,
					sample_id,
					threads
				): chrom for chrom in chromosomes
			}

			for future in as_completed(future_to_chromosome):
				chrom = future_to_chromosome[future]
				try:
					future.result()
				except subprocess.CalledProcessError as e:
					print(f'[ERROR] Failed processing {chrom} for {sample_id}: {e}')
					sys.stdout.flush()
				except Exception as e:
					print(f'[ERROR] Unexpected error for {chrom} in {sample_id}: {e}')
					sys.stdout.flush()

def main():
	parser = argparse.ArgumentParser(
		description='Split BAM files by autosomes in parallel.'
	)
	parser.add_argument('-i', '--input-dir', required=True, help='Path to directory containing BAM files')
	parser.add_argument('-o', '--output-dir', default='.', help='Directory to store output BAM files')
	parser.add_argument('-s', '--samples-file', required=True, help='Path to sample file with genetic IDs')
	parser.add_argument('-t', '--threads', type=int, default=1, help='Threads per samtools call')
	parser.add_argument('-c', '--chromosomes', nargs='+', default=[f'chr{x}' for x in range(1, 23)], help='Chromosomes to process')

	args = parser.parse_args()

	os.makedirs(args.output_dir, exist_ok=True)

	bam_files = get_bam_files(args.input_dir)
	sample_groups = group_samples(args.input_dir, args.samples_file, bam_files)

	for sample_id, bam_files in sample_groups.items():
		print(f'[INFO] Processing sample: {sample_id}')
		sys.stdout.flush()
		process_sample(sample_id, bam_files, args.chromosomes, args.output_dir, args.threads)

	print('[INFO] Splitting complete!')
	sys.stdout.flush()

if __name__ == '__main__':
	main()

