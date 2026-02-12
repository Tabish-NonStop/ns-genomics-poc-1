nextflow.enable.dsl=2

include { FASTQC }                      from './modules/fastqc/main.nf'
include { MULTIQC }                     from './modules/multiqc/main.nf'
include { FASTP_TRIMMING }              from './modules/fastp/trimming/main.nf'
include { BWA_ALIGNMENT }               from './modules/bwa-mem/alignment/main.nf'
include { SAMTOOLS_COMPRESS_AND_SORT }  from './modules/samtools/compress_and_sort/main.nf'
include { SAMTOOLS_INDEXING }           from './modules/samtools/indexing/main.nf'
include { GATK_HAPLOTYPECALLER }        from './modules/gatk/haplotypecaller/variant_calling/main.nf'

include { QUALITY_CONTROL }             from './subworkflows/quality_control/main.nf'
include { VARIANT_CALLING }             from './subworkflows/variant_calling/main.nf'

workflow {

    reads_ch                    = Channel.fromPath(params.reads)

    if (params.run_qc) {
        QUALITY_CONTROL(reads_ch)
    }

    trimmed_reads_ch            = FASTP_TRIMMING(reads_ch)
    trimmed_reads_mapped_ch     = trimmed_reads_ch.map { trimmed, html, json -> trimmed }  

    reference_fasta_ch          = Channel.fromPath(params.reference)
    reference_fai_ch            = Channel.fromPath(params.reference_fai)
    reference_dict_ch           = Channel.fromPath(params.reference_dict)

    if (params.run_vc) {
        _vcf_ch = VARIANT_CALLING(
            reads_ch,
            reference_fasta_ch,
            reference_fai_ch,
            reference_dict_ch
        )
    }
}
