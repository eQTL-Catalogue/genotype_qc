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
    file "output.vcf.gz" into mapped_vcf_ch

    shell:
    """
    CrossMap.py vcf ${chain_file} ${vcf} ${ref_genome} output.vcf
    bgzip output.vcf
    """
}

process filter_vcf{
    input:
    each file(vcf) from mapped_vcf_ch
    val r2_thresh from Channel.from(params.r2_thresh)

    output:
    file "output.vcf.gz" into filtered_vcf_ch

    shell:
    """
    bcftools filter -i 'INFO/R2 > ${r2_thresh}' ${vcf} -Oz -o {output.vcf}
    """
}

