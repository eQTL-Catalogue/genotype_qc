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
    .into { ref_panel_harmonise_genotypes; ref_panel_vcf_fixref } 

process harmonise_genotypes{
    input:
    set file(study_name_bed), file(study_name_bim), file(study_name_fam) from bfile_ch
    set file(vcf_file), file(vcf_file_index) from ref_panel_harmonise_genotypes.collect()

    output:
    set file("harmonised.bed"), file("harmonised.bim"), file("harmonised.fam") into harmonised_genotypes

    script:
    """
    java -jar $baseDir/bin/GenotypeHarmonizer-1.4.20/GenotypeHarmonizer.jar\
     --input ${study_name_bed.baseName}\
     --inputType PLINK_BED\
     --ref ${vcf_file.simpleName}\
     --refType VCF\
     --update-id\
     --output harmonised
    """
}

process plink_to_vcf{

    input:
    set file(bed), file(bim), file(fam) from harmonised_genotypes

    output:
    file "harmonised.vcf.gz" into harmonised_vcf_ch

    script:
    """
    plink2 --bfile ${bed.simpleName} --recode vcf-iid --out ${bed.simpleName}
    bgzip harmonised.vcf
    """
}

process vcf_fixref{
    
    input:
    file input_vcf from harmonised_vcf_ch
    file fasta from ref_genome_ch.collect()
    set file(vcf_file), file(vcf_file_index) from ref_panel_vcf_fixref.collect()

    output:
    file "fixref.vcf.gz" into filter_vcf_input

    script:
    """
    bcftools index ${input_vcf}
    bcftools +fixref ${input_vcf} -- -f ${fasta} -i ${vcf_file} | \
     bcftools norm --check-ref x -f ${fasta} -Oz -o fixref.vcf.gz
    """
}

process filter_vcf{

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: {filename -> if (filename == "filtered.vcf.gz") "${params.output_name}.vcf.gz" else null }

    input:
    file input_vcf from filter_vcf_input

    output:
    set file("filtered.vcf.gz"), file("filtered.vcf.gz.csi") into split_vcf_input, missingness_input

    script:
    """
    #Add tags
    bcftools +fill-tags ${input_vcf} -Oz -o tagged.vcf.gz

    #Filter rare and non-HWE variants and those with abnormal alleles and duplicates
    bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' tagged.vcf.gz |\
     bcftools filter -e 'REF="N" | REF="I" | REF="D"' |\
     bcftools filter -e "ALT='.'" |\
     bcftools norm -d all |\
     bcftools norm -m+any |\
     bcftools view -m2 -M2 -Oz -o filtered.vcf.gz

     #Index the output file
     bcftools index filtered.vcf.gz
    """
}

process calculate_missingness{
    publishDir "${params.outdir}", mode: 'copy',
        saveAs: {filename -> if (filename == "genotypes.imiss") "${params.output_name}.imiss" else null }
    
    input:
    set file(input_vcf), file(input_vcf_index) from missingness_input 

    output:
    file "genotypes.imiss" into missing_individuals

    script:
    """
    vcftools --gzvcf ${input_vcf} --missing-indv --out genotypes
    """
}

process split_by_chr{

    publishDir "${params.outdir}/by_chr", mode: 'copy',
        saveAs: {filename -> if (filename.indexOf(".vcf.gz") > 0) filename else null }
    
    input:
    set file(input_vcf), file(input_vcf_index) from split_vcf_input
    each chr from Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

    output:
    file "chr_${chr}.vcf.gz" into individual_chromosomes

    script:
    """
    bcftools view -r ${chr} ${input_vcf} -Oz -o chr_${chr}.vcf.gz
    """
}