impHet - a pipeline for estimating heterozygosity in aDNA using imputation
===============================================================

Overview
--------
The `impHet` pipeline takes BAM files as input along with a phased reference panel, processes them using GLIMPSE2 to produce imputed genotypes, and computes heterozygosity for genomic regions in the input samples. `impHet` also performs downstream analyses of the heterozygosity estimates it generates: (1) heterozygosity estimates for specific `labeled' regions, such as the MHC, which are added to the input samples metadata file; (2) the results of regressions of heterozygosity over time to test which genomic regions underwent substantial changes in heterozygosity over the entire period encompassed by the input samples; (3) estimation of selection coefficients based on the change in heterozygosity between fixed temporal windows; (4) temporal trends of heterozygosity for genomic regions (intended as a downstream analysis for the results of (2)).

This repository also includes utilities to facilitate use of `impHet`. To prepare the reference, we include a utility to liftover hg38 reference panels to hg19, and generate additional files required by GLIMPSE2. We also include a utility for splitting BAM files into chromosomes, which is a requirement of `impHet`.

Finally, we include the validation pipeline that can be used to compute the expected error of the heterozygosity estimate, given different filtering of depth (for non-imputed genotypes) or genotype probability (for imputed genotypes from GLIMPSE2). We do not provide a minor allele frequency (MAF) filter option, but users can easily perform such filtering downstream of the validation pipeline (using the PLINK files) to assess the effect on the error.


Folder layout
-------------
```
Reference layout
        <reference_root>/<prefix>_<chrom>/
                <prefix>_<chrom>.bcf
                <prefix>_<chrom>.bcf.csi
                <prefix>_<chrom>.sites.vcf.gz
        <reference_root>/map/chr<chrom>.b<N>.gmap.gz (Genetic map, must be downloaded separately)

validation/
    downsampled/
        chr<N>/ SOURCE_TEMPLATE_chr<N>.bam

    called/
        SAMPLE/
            SAMPLE_chr<N>.bcf          per-chromosome BCFs (depth-filtered later)
            SAMPLE.bcf                 concatenated across chromosomes
            SAMPLE_filtered.bcf        +setGT depth-filtered

    truth/
        SAMPLE/
            SAMPLE_chr<N>.bcf
            SAMPLE.bcf
            SAMPLE_filtered.bcf

    imputed/
        SAMPLE/
            SAMPLE_chr<N>.bcf          per-chromosome imputed BCFs
            SAMPLE.bcf
            SAMPLE_filtered.bcf
        merged_full.bcf
        merged_full_filtered.bcf       +setGT GP-filtered

```

Prerequisites
-------------
* python >= 3.11
* bcftools >= 1.17 **compiled with `+liftover` plugin**  
* samtools  
* bedtools  
* plink 1.90 beta
* GLIMPSE2 binaries (`<glimpse_dir>/<stage>/bin/GLIMPSE2_<stage>`)
* Python packages (see `requirements.txt`)

GPU (cupy)
------------
* Note that we currently require the Python package `cupy` (which requires a compatible GPU card) to speedup the bootstrapped regressions. We will write a CPU-based version for this in the near future.

Getting started
----------------
* Clone this repository to your machine
* Install the Python packages using `python3 -m pip install -r requirements.txt`
* Acquire the necessary data and setup input directory structure (see "Input data" below)

Input data
----------
* Phased reference panel (hg19 or hg38) (`*.bcf`)  
* UCSC chain file `hg38ToHg19.over.chain.gz` (only for liftover mode)
* hg38 FASTA and hg19 FASTA
* recombination map (e.g., `map/chr<N>.b37.gmap.gz`)
* sample metadata TSV containing a `date' column

Entry-point scripts
-------------------

1. prepare_reference.py  
   liftover **or** clean one chromosome panel
   ```
   python prepare_reference.py \
       --input-bcf gnomad_chr6.bcf \
       --output-dir reference \
       --contig 6 \
       --mode lift \
       --hg38-fa hg38.fa --hg19-fa hg19.fa \
       --chain-file hg38ToHg19.over.chain.gz
   ```

2. impute_genotypes.py  
   split reference, GLIMPSE2 imputation, PLINK conversion
   ```
   python impute_genotypes.py \
       --source-bam-directory bam_high \
       --reference-directory reference \
       --output-dir results \
       --glimpse-directory /path/to/glimpse2 \
       --threads 16
   ```

3. validation.py
   downsample -> genotype-call -> impute -> merge -> PLINK

   ```
   python validation.py \
       --source-bam-directory   bam_high \
       --template-bam-directory bam_low \
       --reference-prefix       reference/gnomad
       --glimpse-directory      /path/to/glimpse2 \
       --reference-fasta        hg19.fa \
       --output-dir             validation \
       --threads                16 \
       --min-gp                 0.9   \
       --min-dp                 4
   ```

4. compute_heterozygosity.py  
   compute window / gene heterozygosity and run bootstrap regressions
   ```
   python compute_heterozygosity.py \
       --samples-csv samples.tsv \
	--bed-prefix plink/imputed/merged_full \
	--labelled-regions 'mhc=6:29690000-33111000;mhc_class_1=6:29690000-31400000' \
	--regions-file genes.bed \
	--window-size        1000000   \  # genomic window width in bp (for null selection coefficient distribution)
	--time-window-start  12000     \  # first window upper bound (years BP) (for selection coefficient computation)
	--time-window        2000      \  # width of each temporal window (for selection coefficient computation)
	--time-step          1000      \  # slide step between windows (for selection coefficient computation)
	--threads            8         \
	--output-path        heterozygosity

   ```

5. gene_heterozygosity_trends.py  
   heterozygosity trends for significant genes
   ```
   python gene_heterozygosity_trends.py \
       --regressions-tsv heterozygosity/regressions.csv \
       --bed-prefix plink/imputed/merged_full \
       --samples-csv samples.tsv \
       --output-tsv heterozygosity/trends.tsv
   ```

Please note
--------

1. Most aDNA libraries are released in hg19, but a few studies use hg38. `impHet` does not check the assembly of the input BAM files, so check this in advance.
2. Some aDNA libraries are released with non-numeric ("chr") chromosome notation. As with assembly mismatches, `impHet` does not check that chromosome notation matches between the BAMs and reference, so users are advised to check this in advance.
