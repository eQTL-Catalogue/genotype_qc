# Genotype Quality Control
This repository contains three workflows for performing genotype data QC for the eQTL Catalogue project. 

### Dependencies
Most of the software dependencies for the pipelines are listed in the [conda environment](https://github.com/kauralasoo/genotype_qc/blob/master/environment.yml) file. The pipelines also depend on LDAK and [GenotypeHarmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer) and [LDAK5](http://dougspeed.com/ldak/) that need to be downladed separately.

## 1. Pre-imputation QC (pre-imputation.nf)
Preparing genotype data for imputation to the 1000 Genomes Phase 3 reference panel.

QC steps:
- Align raw genotypes to the reference panel with [Genotype Harmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer).
- Convert the genotypes to the VCF format with [PLINK](https://www.cog-genomics.org/plink/1.9/). 
- Exclude variants with Hardy-Weinberg p-value < 1e-6, missingness > 0.05 and minor allele frequency < 0.01 with [bcftools](https://samtools.github.io/bcftools/)
- Calculate individual-level missingness using [vcftools](https://vcftools.github.io/perl_module.html).
- Create separate VCF files for each chromosome.

## 2. Convert imputed genotypes to GRCh38 coordinates (crossmap.nf)

## 3. Project individuals to reference populations from 1000 Genomes Project.

```bash
nextflow run pop_assign.nf -profile pop_assign --vcf <path_to_vcf.vcf.gz> --data_name <study_name>
```


