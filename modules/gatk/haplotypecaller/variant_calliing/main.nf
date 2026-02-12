process GATK_HAPLOTYPECALLER {
    container 'broadinstitute/gatk:latest'

    publishDir "${params.outdir}/gatk/variant_calling/", mode: 'copy'

    input:
    tuple path(sorted_bam), path(sorted_bam_index)
    path reference_fasta
    path reference_fai
    path reference_dict

    output:
    path "${sorted_bam.simpleName}.vcf.gz"

    script:
    """
    gatk HaplotypeCaller \
        -R ${reference_fasta} \
        -I ${sorted_bam} \
        -O ${sorted_bam.simpleName}.vcf.gz \
        --native-pair-hmm-threads 4
    """
}