# Genotype Quality Control
This repository contains three workflows for performing genotype data QC for the eQTL Catalogue project. 

### Dependencies
Most of the software dependencies for the pipelines are listed in the [conda environment](https://github.com/kauralasoo/genotype_qc/blob/master/environment.yml) file. Docker container with all of these dependencies can be obtained from [DockerHub](https://hub.docker.com/r/kauralasoo/genotype_qc).

The pipelines also require [GenotypeHarmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer) and [LDAK5](http://dougspeed.com/ldak/) that need to be downladed separately. Script for downloading those can be found [here](https://github.com/kauralasoo/genotype_qc/blob/master/download_binaries.sh).

## 1. Pre-imputation QC (pre-imputation.nf)
Preparing genotype data for imputation to the 1000 Genomes Phase 3 reference panel.

QC steps:
- Align raw genotypes to the reference panel with [Genotype Harmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer).
- Convert the genotypes to the VCF format with [PLINK](https://www.cog-genomics.org/plink/1.9/). 
- Exclude variants with Hardy-Weinberg p-value < 1e-6, missingness > 0.05 and minor allele frequency < 0.01 with [bcftools](https://samtools.github.io/bcftools/)
- Calculate individual-level missingness using [vcftools](https://vcftools.github.io/perl_module.html).
- Create separate VCF files for each chromosome.

## 2. Convert imputed genotypes to GRCh38 coordinates (crossmap.nf)

## 3. Project individuals to 1000 Genomes Project reference populations (pop_assign.nf).

#### Input

Genotype data imputed to 1000 Genomes Phase 3 reference panel.

#### Analysis steps

- Perform LD pruning on the reference dataset with [PLINK](https://www.cog-genomics.org/plink/1.9/).
- Perform PCA and project new samples to the reference principal components with [LDAK](http://dougspeed.com/ldak/).

```bash
nextflow run pop_assign.nf -profile pop_assign --vcf <path_to_vcf.vcf.gz> --data_name <study_name>
```
Authors: Initial version of the population assignment pipeline was implemented by [Katerina Peikova](https://github.com/peikovakate) and Marija Samoviča, later modified by [Nurlan Kerimov](https://github.com/kerimoff/) and Kaur Alasoo.

