Channel
    .fromPath(params.chain_file)
    .ifEmpty { exit 1, "CrossMap.py chain file not found: ${params.chain_file}" } 
    .set { chain_file_ch }

Channel
    .fromPath(params.ref_genome)
    .ifEmpty { exit 1, "Reference genome fasta file not found: ${params.ref_genome}" } 
    .set { ref_genome_ch }

Channel
    .fromPath( params.vcf_files )
    .ifEmpty { exit 1, "GRCh37 vcf files not found: ${params.vcf_files}" }
    .set { vcf_file_ch }

process crossmap_genotypes{
    input:
    file chain_file from chain_file_ch
    file ref_genome from ref_genome_ch
    each file(vcf) from vcf_file_ch

    output:
    file "${vcf.simpleName}_mapped.vcf.gz" into mapped_vcf_ch

    shell:
    """
    CrossMap.py vcf ${chain_file} ${vcf} ${ref_genome} ${vcf.simpleName}_mapped.vcf
    bgzip ${vcf.simpleName}_mapped.vcf
    """
}

process filter_vcf{
    input:
    each file(vcf) from mapped_vcf_ch
    val r2_thresh from Channel.from(params.r2_thresh)

    output:
    file "${vcf.simpleName}.vcf.gz" into filtered_vcf_ch

    shell:
    """
    bcftools filter -i 'INFO/R2 > ${r2_thresh}' ${vcf} -Oz -o ${vcf.simpleName}.vcf.gz
    """
}

process merge_vcf{
    input:
    file input_files from filtered_vcf_ch.collect()

    output:
    file "output.vcf.gz" into merged_vcf_ch

    shell:
    """
    bcftools concat ${input_files} | bcftools sort -Oz -o output.vcf.gz
    """
}

process keep_chromosomes{
    input:
    file input_vcf from merged_vcf_ch

    output:
    file "output.vcf.gz" into final_vcf_ch

    shell:
    """
    bcftools index ${input_vcf}
    bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X ${input_vcf} |\
     bcftools annotate --set-id 'chr%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -Oz -o output.vcf.gz
    """
}

