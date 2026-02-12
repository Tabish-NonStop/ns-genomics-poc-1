#! usr/bin/env nextflow
nextflow.enable.dsl=2

process SAMTOOLS_INDEXING {
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    
    publishDir "${params.outdir}/indexing/samtools", mode: 'copy'

    input: 
    path bam_file

    output: 
    tuple path (bam_file), path ("${bam_file}.bai")

    script:
    """
    samtools index ${bam_file}
    """
}