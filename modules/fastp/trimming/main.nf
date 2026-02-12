#! usr/bin/env nextflow
nextflow.enable.dsl=2

process FASTP_TRIMMING {
    container 'biocontainers/fastp:v0.20.1_cv1'

    publishDir "${params.outdir}/trimming/fastp", mode: 'copy'

    input:
    path reads

    output:
    tuple path("${reads.simpleName}.trimmed.fastq.gz"),
          path("${reads.simpleName}.fastp.html"), 
          path("${reads.simpleName}.fastp.json")

    script:
    """
    fastp \
        --in1 ${reads} \
        --out1 ${reads.simpleName}.trimmed.fastq.gz \
        --html ${reads.simpleName}.fastp.html \
        --json ${reads.simpleName}.fastp.json \
        --qualified_quality_phred 20 \
        --length_required 50
    """
}
