'''
GLIMPSE2 Python wrapper

Public entry points
	prepare_reference(glimpse_dir, reference_prefix_path, chromosome)
	impute_sample(glimpse_dir, reference_prefix_path, chromosome, bam_file, work_root)

Reference layout
	<reference_root>/<prefix>_<chrom>/
		<prefix>_<chrom>.bcf
		<prefix>_<chrom>.bcf.csi
		<prefix>_<chrom>.sites.vcf.gz

Genetic map
	<reference_root>/map/chr<chrom>.b37.gmap.gz
'''

import logging
import os
import re
import subprocess
import re as _re_nat
from pathlib import Path


def _init_logger():
	log_path = Path.cwd() / 'glimpse2_wrapper.log'
	log_path.parent.mkdir(parents=True, exist_ok=True)
	logger = logging.getLogger('glimpse2.wrapper')
	if not logger.handlers:
		logger.setLevel(logging.INFO)
		fmt = logging.Formatter('[%(levelname)s] [GLIMPSE2 -> %(stage)s -> %(sample)s]: %(message)s')
		handler = logging.FileHandler(log_path)
		handler.setFormatter(fmt)
		logger.addHandler(handler)

		stdout_path = Path.cwd() / 'glimpse2_stdout.log'
		stdout_handler = logging.FileHandler(stdout_path)
		stdout_handler.setFormatter(logging.Formatter('[%(stage)s] [%(sample)s]\n%(message)s'))
		stdout_handler.setLevel(logging.DEBUG)
		stdout_handler.addFilter(lambda rec: rec.levelno == logging.DEBUG)
		logger.addHandler(stdout_handler)
	return logger


def _run(cmd, logger, stage, sample):
	process = subprocess.run(cmd, text=True, capture_output=True)
	logger.debug(process.stdout, extra={'stage': stage, 'sample': sample})
	if process.stderr:
		logger.info('stderr:\n%s', process.stderr.strip(), extra={'stage': stage, 'sample': sample})
	if process.returncode:
		raise RuntimeError(f'{stage} failed for {sample}.')


def _validate_chromosome(label):
	if not re.fullmatch(r'([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)', label):
		raise ValueError(f'Unsupported chromosome label: {label}')


def _index_exists(bcf_path):
	for suff in ('.csi', '.tbi'):
		if bcf_path.with_suffix(bcf_path.suffix + suff).exists():
			return True
	return False


class _GlimpseBase:
	def __init__(self, glimpse_dir, reference_prefix_path, chromosome, logger=None):
		_validate_chromosome(chromosome)
		self.gl_dir = Path(glimpse_dir).resolve()
		self.prefix_path = Path(reference_prefix_path).resolve()
		self.prefix_name = self.prefix_path.name
		self.ref_root = self.prefix_path.parent
		self.chrom = chromosome

		self.ref_chr_dir = self.ref_root / f'{self.prefix_name}_{self.chrom}'
		self.bcf = self.ref_chr_dir / f'{self.prefix_name}_{self.chrom}.bcf'
		self.bcf_sites = self.ref_chr_dir / f'{self.prefix_name}_{self.chrom}.sites.vcf.gz'
		self.chunks_txt = self.ref_chr_dir / f'{self.prefix_name}_{self.chrom}.chunks.txt'
		self.map_file = self.ref_root / 'map' / f'chr{self.chrom}.b37.gmap.gz'

		if not self.bcf.exists() or not self.bcf_sites.exists():
			raise FileNotFoundError(f'Reference files missing in {self.ref_chr_dir}')
		if not _index_exists(self.bcf):
			raise FileNotFoundError(f'BCF index missing for {self.bcf}')

		self.logger = logger or _init_logger()

	def exe(self, sub_path):
		return self.gl_dir / sub_path


class ReferencePreparer(_GlimpseBase):
	def run_chunk(self):
		if self.chunks_txt.exists():
			self.logger.info('%s already present – skipping chunk step', self.chunks_txt, extra={'stage': 'chunk', 'sample': 'reference'})
			return
		cmd = [
			str(self.exe('chunk/bin/GLIMPSE2_chunk')),
			'--input', str(self.bcf_sites),
			'--region', self.chrom,
			'--output', str(self.chunks_txt),
			'--map', str(self.map_file),
			'--sequential',
		]
		_run(cmd, self.logger, 'chunk', 'reference')
		if not self.chunks_txt.exists():
			raise RuntimeError('CHUNK produced no output.')

	def run_split_reference(self):
		split_dir = self.ref_chr_dir / 'split'
		split_dir.mkdir(exist_ok=True)
		exe_path = self.exe('split_reference/bin/GLIMPSE2_split_reference')
		ref_prefix = split_dir / f'{self.prefix_name}'

		with self.chunks_txt.open() as fh:
			for line in fh:
				_, _, irg, org, *_ = line.split()
				start, end = irg.split(':')[1].split('-')
				bin_path = ref_prefix.with_name(f'{self.prefix_name}_{self.chrom}_{start}_{end}.bin')
				if bin_path.exists():
					continue
				cmd = [
					str(exe_path),
					'--reference', str(self.bcf),
					'--map', str(self.map_file),
					'--input-region', irg,
					'--output-region', org,
					'--output', str(ref_prefix),
				]
				_run(cmd, self.logger, 'split_reference', 'reference')
				if not bin_path.exists():
					raise RuntimeError('SPLIT_REFERENCE failed to create ' + str(bin_path))

	def prepare(self):
		self.run_chunk()
		self.run_split_reference()


class SampleImputer(_GlimpseBase):
	def __init__(self, glimpse_dir, reference_prefix_path, chromosome, bam_file, work_root=None, logger=None):
		super().__init__(glimpse_dir, reference_prefix_path, chromosome, logger)
		self.bam = Path(bam_file).resolve()
		if not self.bam.exists():
			raise FileNotFoundError(self.bam)
		self.sample_id = self.bam.stem
		work_root = Path(work_root) if work_root else Path(os.getenv('TMPDIR', '/tmp'))
		self.work_dir = work_root.resolve() / self.sample_id
		self.phase_dir = self.work_dir / f'GLIMPSE_impute_chr{self.chrom}'
		self.work_dir.mkdir(parents=True, exist_ok=True)
		self.phase_dir.mkdir(parents=True, exist_ok=True)
		self.split_dir = self.ref_chr_dir / 'split'
		if not self.split_dir.exists():
			raise FileNotFoundError('Split reference not found – run prepare_reference first.')

	def run_phase(self):
		exe_phase = self.exe('phase/bin/GLIMPSE2_phase')
		for bin_file in sorted(self.split_dir.glob('*.bin')):
			out_bcf = self.phase_dir / f'imputed_{bin_file.stem}.bcf'
			if out_bcf.exists():
				continue
			cmd = [
				str(exe_phase),
				'--bam-file', str(self.bam),
				'--reference', str(bin_file),
				'--output', str(out_bcf),
			]
			_run(cmd, self.logger, 'phase', self.sample_id)
			if not out_bcf.exists():
				raise RuntimeError('PHASE failed to create ' + str(out_bcf))

	def run_ligate(self):
		def _nat_key(p):
			return [int(tok) if tok.isdigit() else tok for tok in _re_nat.split(r'(\d+)', p.stem)]

		lst_path = self.work_dir / f'list_chr{self.chrom}.txt'
		with lst_path.open('w') as fh:
			for bcf_file in sorted(self.phase_dir.glob('imputed_*.bcf'), key=_nat_key):
				fh.write(str(bcf_file) + '\n')

		ligated_bcf = self.work_dir / f'ligated_chr{self.chrom}.bcf'
		if ligated_bcf.exists():
			return
		cmd = [
			str(self.exe('ligate/bin/GLIMPSE2_ligate')),
			'--input', str(lst_path),
			'--output', str(ligated_bcf),
		]
		_run(cmd, self.logger, 'ligate', self.sample_id)
		if not ligated_bcf.exists():
			raise RuntimeError('LIGATE failed for ' + self.sample_id)

	def impute(self):
		self.run_phase()
		self.run_ligate()
		return self.work_dir / f'ligated_chr{self.chrom}.bcf'


def prepare_reference(glimpse_dir, reference_prefix_path, chromosome):
	ReferencePreparer(glimpse_dir, reference_prefix_path, chromosome).prepare()


def impute_sample(glimpse_dir, reference_prefix_path, chromosome, bam_file, work_root=None):
	return SampleImputer(glimpse_dir, reference_prefix_path, chromosome, bam_file, work_root).impute()
