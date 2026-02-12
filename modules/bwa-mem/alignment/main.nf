#! usr/bin/env nextflow
nextflow.enable.dsl=2

process BWA_ALIGNMENT {
    container 'biocontainers/bwa:v0.7.17_cv1'

    publishDir "${params.outdir}/alignment/bwa", mode: 'copy'

    input:
    path trimmed_reads
    path reference_fasta

    output:
    path "${trimmed_reads.baseName}.sam"

   script:
    """
    # Define bash variable
    SAMPLE=${trimmed_reads.baseName}

    # Index reference (will be cached by -resume)
    bwa index ${reference_fasta}

    bwa mem -t 4 \
        -R "@RG\\tID:\${SAMPLE}\\tSM:\${SAMPLE}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1" \
        ${reference_fasta} \
        ${trimmed_reads} \
        > \${SAMPLE}.sam
    """
}