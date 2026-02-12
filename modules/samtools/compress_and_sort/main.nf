#! usr/bin/env nextflow
nextflow.enable.dsl=2

process SAMTOOLS_COMPRESS_AND_SORT {
    container 'biocontainers/samtools:v1.9-4-deb_cv1'

    publishDir "${params.outdir}/compression-and-sorting/samtools", mode: 'copy'

    input:
    path sam_file

    output:
    path "${sam_file.baseName}.sorted.bam"

    script:
    """
    samtools sort -o ${sam_file.baseName}.sorted.bam ${sam_file}
    """
}