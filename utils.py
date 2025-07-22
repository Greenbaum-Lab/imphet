import subprocess

def create_bcf_from_bam(input_bam, chromosome, reference_fasta, reference_sites, output_path, tsv=None):
	'''
	This function creates a BCF file from a given BAM file, using the same approach that appears in the GLIMPSE2 manual for testing concordance
	It runs mpileup and then calls bcftools index on the resulting output.

	Parameters:
	- input_bam: Path to the input BAM file.
	- chromosome: The chromosome for mpileup.
	- reference_fasta: The reference FASTA file for base calling.
	- reference_sites: A file of reference sites for variant calling.
	- output_path: Path to the output BCF file.
	- tsv: Optional TSV file containing allele annotations.
	'''
	try:
		command = subprocess.run(f'bcftools mpileup '
			f'-f {reference_fasta} '
			f'--ignore-RG '	# Ignore read groups in BAMs, this mirrors GLIMPSE2 behavior with respect to read groups
			f'-I '
			f'-E '
			f'-a 'FORMAT/DP' '
			f'-T {reference_sites} ' # Call variants from reference panel only
			f'-r {chromosome} '
			f'{input_bam} -Ou | '
			f'bcftools call -Aim '
			f'-C alleles '
			f'-T {tsv} ' # Call variants from reference panel only
			f'-Ob -o {output_path}',
			shell=True, capture_output=True, text=True)

		if command.returncode:
			raise Exception(f'bcftools mpileup failed: {command.stderr}')

		command = subprocess.run(['bcftools', 'index', output_path], capture_output=True, text=True)
		if command.returncode:
			raise Exception(f'bcftools index failed: {command.stderr}')

	except Exception as e:
		print(f'Error: {str(e)}')
		raise
