#!/usr/bin/env python3

'''
Prepare a chromosome specific reference panel.

Modes
	lift  : liftover hg38 panel to hg19 then clean
	clean : clean an existing hg19 panel

Cleaning steps
	1. keep bi-allelic SNPs
	2. remove duplicate positions
	3. left-align and normalise against hg19 fasta
	4. rename chromosomes (chr1 → 1, chrM → MT)
	5. keep requested contig only
	6. index final BCF
	7. produce .sites.vcf.gz and .sites.tsv.gz with indexes

Output layout
	<output_dir>/<prefix>_<contig>/
		<prefix>_<contig>.bcf
		<prefix>_<contig>.bcf.csi
		<prefix>_<contig>.sites.vcf.gz
		<prefix>_<contig>.sites.vcf.gz.csi
		<prefix>_<contig>.sites.tsv.gz
		<prefix>_<contig>.sites.tsv.gz.tbi
'''

import argparse
import subprocess
import tempfile
import shlex
from datetime import datetime
from pathlib import Path


def _run_cmd(cmd, log_path):
	with open(log_path, 'a') as log_fh:
		log_fh.write(f'[{datetime.now().isoformat()}] {cmd}\n')
		subprocess.check_call(cmd, shell=True, stderr=log_fh, stdout=log_fh)


def _index_bcf(bcf_path, threads, log_path):
	idx_path = bcf_path.with_suffix(bcf_path.suffix + '.csi')
	if idx_path.exists():
		return
	_run_cmd(f'bcftools index -f --threads {threads} {shlex.quote(str(bcf_path))}', log_path)


def _filter_snps(in_bcf, out_bcf, threads, log_path):
	if out_bcf.exists():
		return
	_run_cmd(
		f'bcftools view -m2 -M2 -v snps --threads {threads} -Ob '
		f'-o {shlex.quote(str(out_bcf))} {shlex.quote(str(in_bcf))}',
		log_path
	)


def _deduplicate(in_bcf, out_bcf, threads, log_path):
	if out_bcf.exists():
		return
	_run_cmd(
		f'bcftools norm -d both --threads {threads} -Ob '
		f'-o {shlex.quote(str(out_bcf))} {shlex.quote(str(in_bcf))}',
		log_path
	)


def _liftover(in_bcf, chain_file, src_fa, dst_fa, out_bcf, log_path):
	if out_bcf.exists():
		return
	sort_tmp = tempfile.mkdtemp(prefix='bcsort_')
	cmd = (
		f'bcftools +liftover --no-version -Ou {shlex.quote(str(in_bcf))} -- '
		f'-s {shlex.quote(str(src_fa))} -f {shlex.quote(str(dst_fa))} '
		f'-c {shlex.quote(str(chain_file))} --write-src | '
		f'bcftools sort -T {sort_tmp} -Ob -o {shlex.quote(str(out_bcf))}'
	)
	_run_cmd(cmd, log_path)


def _left_align(in_bcf, ref_fa, out_bcf, threads, log_path):
	if out_bcf.exists():
		return
	_run_cmd(
		f'bcftools norm -f {shlex.quote(str(ref_fa))} -d both --threads {threads} '
		f'-Ob -o {shlex.quote(str(out_bcf))} {shlex.quote(str(in_bcf))}',
		log_path
	)


def _rename_chroms(in_bcf, out_bcf, threads, log_path):
	if out_bcf.exists():
		return
	with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
		map_path = Path(tmp.name)
		for i in range(1, 23):
			tmp.write(f'chr{i}\t{i}\n')
		tmp.write('chrX\tX\nchrY\tY\nchrM\tMT\n')
	_run_cmd(
		f'bcftools annotate --rename-chrs {map_path} --threads {threads} '
		f'-Ob -o {shlex.quote(str(out_bcf))} {shlex.quote(str(in_bcf))}',
		log_path
	)
	map_path.unlink(missing_ok=True)


def _restrict_contig(in_bcf, contig, out_bcf, threads, log_path):
	if out_bcf.exists():
		return
	_run_cmd(
		f'bcftools view -r {shlex.quote(contig)} --threads {threads} -Ob '
		f'-o {shlex.quote(str(out_bcf))} {shlex.quote(str(in_bcf))}',
		log_path
	)


def _create_site_indexes(final_bcf, threads, log_path):
	site_vcf = final_bcf.with_suffix('.sites.vcf.gz')
	if site_vcf.exists():
		return
	_run_cmd(
		f'bcftools view -G -m2 -M2 -v snps --threads {threads} {shlex.quote(str(final_bcf))} '
		f'-Oz -o {shlex.quote(str(site_vcf))}',
		log_path
	)
	_index_bcf(site_vcf, threads, log_path)
	site_tsv = final_bcf.with_suffix('.sites.tsv.gz')
	_run_cmd(
		f"bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' {site_vcf} | "
		f'bgzip -c > {shlex.quote(str(site_tsv))}',
		log_path
	)
	_run_cmd(f'tabix -s1 -b2 -e2 {shlex.quote(str(site_tsv))}', log_path)


def prepare_reference(args):
	log_path = args.output_prefix.with_suffix('.prep.log')

	snps_bcf = args.output_prefix.with_suffix('.snps.bcf')
	_filter_snps(args.input_bcf, snps_bcf, args.threads, log_path)

	dedup_bcf = args.output_prefix.with_suffix('.dedup.bcf')
	_deduplicate(snps_bcf, dedup_bcf, args.threads, log_path)

	step_bcf = dedup_bcf
	if args.mode == 'lift':
		lift_bcf = args.output_prefix.with_suffix('.lift.bcf')
		_liftover(dedup_bcf, args.chain_file, args.hg38_fa, args.hg19_fa, lift_bcf, log_path)
		step_bcf = lift_bcf

	norm_bcf = args.output_prefix.with_suffix('.norm.bcf')
	_left_align(step_bcf, args.hg19_fa, norm_bcf, args.threads, log_path)

	rename_bcf = args.output_prefix.with_suffix('.renamed.bcf')
	_rename_chroms(norm_bcf, rename_bcf, args.threads, log_path)
	_index_bcf(rename_bcf, args.threads, log_path)

	final_bcf = args.output_prefix.with_suffix('.bcf')
	_restrict_contig(rename_bcf, args.contig, final_bcf, args.threads, log_path)
	_index_bcf(final_bcf, args.threads, log_path)

	_create_site_indexes(final_bcf, args.threads, log_path)

	for tmp_path in [snps_bcf, dedup_bcf, norm_bcf, rename_bcf]:
		tmp_path.unlink(missing_ok=True)
		Path(str(tmp_path) + '.csi').unlink(missing_ok=True)
	if args.mode == 'lift':
		lift_tmp = args.output_prefix.with_suffix('.lift.bcf')
		lift_tmp.unlink(missing_ok=True)
		Path(str(lift_tmp) + '.csi').unlink(missing_ok=True)


def main():
	parser = argparse.ArgumentParser(description='Prepare reference panel chromosome')
	parser.add_argument('--input-bcf', required=True)
	parser.add_argument('--output-dir', required=True)
	parser.add_argument('--prefix', required=True)
	parser.add_argument('--contig', required=True)
	parser.add_argument('--hg38-fa')
	parser.add_argument('--hg19-fa', required=True)
	parser.add_argument('--chain-file')
	parser.add_argument('--threads', type=int, default=4)
	parser.add_argument('--mode', choices=['lift', 'clean'], default='clean')
	args = parser.parse_args()

	sub_dir = Path(args.output_dir) / f'{args.prefix}_{args.contig}'
	sub_dir.mkdir(parents=True, exist_ok=True)
	args.output_prefix = sub_dir / f'{args.prefix}_{args.contig}'

	if args.mode == 'lift' and (not args.chain_file or not args.hg38_fa):
		raise ValueError('lift mode requires --chain-file and --hg38-fa')

	prepare_reference(args)


if __name__ == '__main__':
	main()
