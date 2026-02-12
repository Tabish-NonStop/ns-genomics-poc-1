#! usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads = "${projectDir}/data/sample_data.fastq"
params.outdir = "${projectDir}/results"
params.reference = "${projectDir}/data/ref/ref.fasta"
params.reference_fai = "${projectDir}/data/ref/ref.fasta.fai"
params.reference_dict = "${projectDir}/data/ref/ref.dict"

include { FASTQC }                      from './modules/fastqc/main.nf'
include { MULTIQC }                     from './modules/multiqc/main.nf'
include { FASTP_TRIMMING }              from './modules/fastp/trimming/main.nf'
include { BWA_ALIGNMENT }               from './modules/bwa-mem/alignment/main.nf'
include { SAMTOOLS_COMPRESS_AND_SORT }  from './modules/samtools/compress_and_sort/main.nf'
include { SAMTOOLS_INDEXING }           from './modules/samtools/indexing/main.nf'





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


workflow {

    reads_ch                    = Channel.fromPath(params.reads)

    fastqc_reports_ch           = FASTQC(reads_ch)

    MULTIQC(fastqc_reports_ch.collect())

    trimmed_reads_ch            = FASTP_TRIMMING(reads_ch)
    trimmed_reads_mapped_ch     = trimmed_reads_ch.map { trimmed, html, json -> trimmed }  

    reference_fasta_ch          = Channel.fromPath(params.reference)
    reference_fai_ch            = Channel.fromPath(params.reference_fai)
    reference_dict_ch           = Channel.fromPath(params.reference_dict)

    aligned_bam_ch              = BWA_ALIGNMENT(trimmed_reads_mapped_ch, reference_fasta_ch)

    compressed_and_sorted_ch    = SAMTOOLS_COMPRESS_AND_SORT(aligned_bam_ch)

    indexed_bam_ch              = SAMTOOLS_INDEXING(compressed_and_sorted_ch)

    vcf_ch = GATK_HAPLOTYPECALLER(
        indexed_bam_ch,
        reference_fasta_ch,
        reference_fai_ch,
        reference_dict_ch
    )

}