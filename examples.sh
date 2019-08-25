#Prepare raw genotype data for imputation
nextflow run pre-imputation_qc.nf -profile eqtl_catalogue -resume\
 --bfile /gpfs/hpc/home/a72094/datasets/controlled_access/Schmiedel_2018/genotypes/matrix/IC_DNA\
 --output_name Schmiedel_2018_GRCh37_genotyped\
 --outdir Schmiedel_2018

#CrossMap imputed VCF files to new genome version
nextflow run crossmap_genotypes.nf -profile crossmap -resume\
 --vcf_files "/gpfs/hpc/home/a72094/datasets/controlled_access/Schmiedel_2018/genotypes/Michigan_GRCh37_Phase3_210819/GRCh37/chr*.dose.vcf.gz"\
 --output_name Schmiedel_2018
