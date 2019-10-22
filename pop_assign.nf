#!/usr/bin/env nextflow

vcf_file = file(params.vcf)
main_vcf = file(params.main_vcf)
populations_file = file(params.populations_file)

// TODO:    
// manage arguments and used programs
// add qtltools option
// add option for only pca without mapping

process main_vcf_to_binary{
    storeDir "$baseDir/reference_vcf_cache/"
    
    input:
    file "vcf_main.vcf.gz" from main_vcf

    output:
    set file('main.bed'), file('main.bim'), file('main.fam') into main_to_delete_dublicates, calculate_relatedness_ch

    script:
    """
    # do ld pruning and save resulted snps in file  
    plink2 --vcf vcf_main.vcf.gz --vcf-half-call h --indep-pairwise 50000 200 0.05 --out main_pruned_varaints_list --threads ${task.cpus}

    # make bfiles for pruned 1000 genome proj 
    plink2 --vcf vcf_main.vcf.gz --vcf-half-call h --extract main_pruned_varaints_list.prune.in --make-bed --out main
    """
}

process calculate_relatedness_matrix{
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set file('main.bed'), file('main.bim'), file('main.fam') from calculate_relatedness_ch
    
    output:
    file("relatedness_matrix.tsv") into relatedness_matrix_ch

    script:
    """
    #Calculate relatedness
    plink2 --make-rel square --bfile main

    #Format relatedness matrix
    Rscript $baseDir/bin/format_kinship.R \\
        --kinship main.rel \\
        --fam main.fam \\
        --out relatedness_matrix.tsv
    """
}

if(params.exclude_population){
    process remove_family{
    storeDir "$baseDir/reference_vcf_cache/"

    input:
    set file('main.bed'), file('main.bim'), file('main.fam') from main_to_delete_dublicates
    file 'ids_to_remove.txt' from file(params.ids_to_remove_file)

    output:
    set file('main_no_dubl.bed'), file('main_no_dubl.bim'), file('main_no_dubl.fam') into  main_to_extract_snps, main_binary_source

    script:
    """
    plink2 --bfile main --remove-fam ids_to_remove.txt --make-bed --out main
    
    # finds dublicate vars
    plink2 --bfile main --list-duplicate-vars --out dubl

    # delete dublicate vars
    plink2 --bfile main --exclude dubl.dupvar --snps-only --make-bed --out main_no_dubl
    """
    }
}else{
    process remove_dubl{
    storeDir "$baseDir/reference_vcf_cache/"

    input:
    set file('main.bed'), file('main.bim'), file('main.fam') from main_to_delete_dublicates
  
    output:
    set file('main_no_dubl.bed'), file('main_no_dubl.bim'), file('main_no_dubl.fam') into main_to_extract_snps, main_binary_source

    script:
    """
    # finds dublicate vars
    plink2 --bfile main --list-duplicate-vars --out dubl
    
    # delete dublicate vars
    plink2 --bfile main --exclude dubl.dupvar --snps-only --make-bed --out main_no_dubl
    """
    }

}

process get_SNPs_list_from_main_dataset{
    storeDir "$baseDir/reference_vcf_cache/"

    input: 
    set file('main_source.bed'), file('main_source.bim'), file('main_source.fam') from main_to_extract_snps

    output:
    file 'main_snps_list.snplist' into main_snps_file

    script:
    """
    plink2 --bfile main_source --write-snplist --out main_snps_list --snps-only
    """
}

// convert vcf file to plink binary file (.bed)
process convertVCFtoBED{
    input:
    file 'source.vcf.gz' from vcf_file

    output:
    set file ('binary_source.bed'), file('binary_source.bim'), file('binary_source.fam') into bed_files

    script:
    """
    plink2 --vcf source.vcf.gz --out binary_source --threads ${task.cpus}
    plink2 --bfile binary_source --list-duplicate-vars --out list_dubl
    plink2 --bfile binary_source --exclude list_dubl.dupvar --snps-only --make-bed --out binary_source
    """
}

process extract_SNPs_and_make_bed{
    input:
    set file ('binary_source.bed'), file('binary_source.bim'), file('binary_source.fam') from bed_files
    file 'main.snplist' from main_snps_file

    output:
    set file('new_dataset_overlapped.bed'), file('new_dataset_overlapped.bim'), file('new_dataset_overlapped.fam') into new_dataset_overlapped1, new_dataset_overlapped2
    file 'overlapped_snps.snplist' into new_dataset_overlapped_snplist

    // extract snps present in pruned data from new dataset
    // --make-bed makes sure bfiles created!
    script:
    """
    plink2 --bfile binary_source --extract main.snplist --make-bed --out new_dataset_overlapped
    plink2 --bfile new_dataset_overlapped --write-snplist --out overlapped_snps
    """
}

process extract_overlapped_SNPs_from_main_dataset {
    input:
    set file('new_dataset_overlapped.bed'), file('new_dataset_overlapped.bim'), file('new_dataset_overlapped.fam') from new_dataset_overlapped1
    file 'overlapped_snps.snplist' from new_dataset_overlapped_snplist
    set file('main_source.bed'), file('main_source.bim'), file('main_source.fam') from main_binary_source

    output:
    set file('main_overlapped.bed'), file('main_overlapped.bim'), file('main_overlapped.fam') into main_to_kins
    set file('main_overlapped.bed'), file('main_overlapped.bim'), file('main_overlapped.fam') into main_bed_to_pca

    script:
    """
    plink2 --bfile main_source --extract overlapped_snps.snplist --make-bed --out main_overlapped
    """
}

process calc_kins_matrices{
    input:
    set file('main_overlapped.bed'), file('main_overlapped.bim'), file('main_overlapped.fam') from main_to_kins

    output:
    set file("main_overlapped_kins.grm.bin"), file("main_overlapped_kins.grm.id"), file("main_overlapped_kins.grm.adjust"), file("main_overlapped_kins.grm.details") into main_to_pca
    
    script:
    """
    ldak --calc-kins-direct main_overlapped_kins --bfile main_overlapped --ignore-weights YES --power -0.25
    """
}

process calc_pca_and_loads{
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set file("main_overlapped_kins.grm.bin"), file("main_overlapped_kins.grm.id"), file("main_overlapped_kins.grm.adjust"), file("main_overlapped_kins.grm.details") from main_to_pca
    set file('main_overlapped.bed'), file('main_overlapped.bim'), file('main_overlapped.fam') from main_bed_to_pca

    output:
    file 'main_overlapped_loads.load' into loads_for_mapping
    file 'main_overlapped_pca.vect' into plot_data_vect
     
    script: 
    """
    ldak --pca main_overlapped_pca --grm main_overlapped_kins --axes $params.num_pc
    ldak --calc-pca-loads main_overlapped_loads --grm main_overlapped_kins --pcastem main_overlapped_pca --bfile main_overlapped
    """
}

process map_new_dataset{
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    set file('new_dataset_overlapped.bed'), file('new_dataset_overlapped.bim'), file('new_dataset_overlapped.fam') from new_dataset_overlapped2
    file 'main_overlapped_loads.load' from loads_for_mapping

    output:
    file 'new_dataset_scores.profile.adj' into plot_data

    script:
    """
    ldak  --calc-scores new_dataset_scores --bfile new_dataset_overlapped --scorefile main_overlapped_loads.load --power 0
    """

}

process plot_pca{
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    file 'new_dataset_scores.profile.adj' from plot_data
    file 'main_overlapped_pca.vect' from plot_data_vect
    file 'samples_data.tsv' from populations_file

    output:
    set file('main_pca.png'), file('projections_only.png'), file('projections_on_main.png'), file('populations.tsv'), file('knn_threshold.png'), file('knn.png')

    script:
    """
    Rscript $baseDir/bin/plot_pca.R main_overlapped_pca.vect new_dataset_scores.profile.adj samples_data.tsv $params.data_name
    """

}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops ... something went wrong" )
}






