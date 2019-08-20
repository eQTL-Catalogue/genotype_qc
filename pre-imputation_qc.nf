Channel
    .fromPath(params.ref_genome)
    .ifEmpty { exit 1, "Reference genome fasta file not found: ${params.ref_genome}" } 
    .set { ref_genome_ch }

Channel
    .from(params.bfile)
    .map { study -> [file("${study}.bed"), file("${study}.bim"), file("${study}.fam")]}
    .set { bfile_ch }

Channel
    .from(params.ref_panel)
    .map { ref -> [file("${ref}.vcf.gz"), file("${ref}.vcf.gz.tbi")]}
    .set { ref_panel_ch } 

process harmonise_genotypes{
    input:
    set file(study_name_bed), file(study_name_bim), file(study_name_fam) from bfile_ch
    set file(vcf_file), file(vcf_file_index) from ref_panel_ch

    output:
    set file("${study_name_bed.simpleName}_harmonised.bed"), file("${study_name_bed.simpleName}_harmonised.bim"), file("${study_name_bed.simpleName}_harmonised.fam") into harmonised_genotypes

    script:
    """
    java -jar $baseDir/bin/GenotypeHarmonizer-1.4.20/GenotypeHarmonizer.jar\
     --input ${study_name_bed.simpleName}\
     --inputType PLINK_BED\
     --ref ${vcf_file.simpleName}\
     --refType VCF\
     --update-id\
     --output ${study_name_bed.simpleName}_harmonised
    """
}

process plink_to_vcf{

    input:
    set file(bed), file(bim), file(fam) from harmonised_genotypes

    output:
    file("${bed.simpleName}.vcf.gz") into harmonised_vcf_ch

    script:
    """
    plink --bfile ${bed.simpleName} --recode vcf-iid --out ${bed.simpleName}
    bgzip ${bed.simpleName}.vcf
    """
}