Channel
    .fromPath(params.ref_genome)
    .ifEmpty { exit 1, "Reference genome fasta file not found: ${params.ref_genome}" } 
    .set { ref_genome_ch }
Channel
    .fromPath(params.dbsnp_vcf)
    .ifEmpty { exit 1, "dbSNP reference vcf file not found: ${params.dbsnp_vcf}" }
    .set { dbsnp_vcf_ch }

Channel
    .from(params.plink_prefix)
    .map { study -> [file("${study}.bed"), file("${study}.bim"), file("${study}.fam")]}
    .set { plink_data }


process harmonise_genotypes{
    input:
    set file(study_name_bed), file(study_name_bim), file(study_name_fam) from plink_data

    output:
    set file("${study_name_bed.simpleName}.harmonised.bed"), file("${study_name_bed.simpleName}.harmonised.bim"), file("${study_name_bed.simpleName}.harmonised.fam") into harmonised_genotypes

    script:
    """
    java -jar $baseDir/bin/GenotypeHarmonizer-1.4.20/GenotypeHarmonizer.jar\
     --input ${study_name_bed.simpleName}\
     --inputType PLINK_BED\
     --ref /gpfs/hpc/home/a72094/datasets/1000G/GRCh37_allele_frequencies\
     --refType VCF\
     --update-id\
     --output ${study_name_bed.simpleName}.harmonized
    """
}

