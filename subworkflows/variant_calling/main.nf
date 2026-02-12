nextflow.enable.dsl = 2

include { FASTP_TRIMMING }              from '../../modules/fastp/trimming/main.nf'
include { BWA_ALIGNMENT }               from '../../modules/bwa-mem/alignment/main.nf'
include { SAMTOOLS_COMPRESS_AND_SORT }  from '../../modules/samtools/compress_and_sort/main.nf'
include { SAMTOOLS_INDEXING }           from '../../modules/samtools/indexing/main.nf'
include { GATK_HAPLOTYPECALLER }        from '../../modules/gatk/haplotypecaller/variant_calling/main.nf'

workflow VARIANT_CALLING {

    take:
        reads_ch
        reference_fasta_ch
        reference_fai_ch
        reference_dict_ch

    main:
    
        trimmed_reads_ch            = FASTP_TRIMMING(reads_ch)
        trimmed_reads_mapped_ch     = trimmed_reads_ch.map { trimmed, _html, _json -> trimmed }  
        aligned_bam_ch              = BWA_ALIGNMENT(trimmed_reads_mapped_ch, reference_fasta_ch)
        compressed_and_sorted_ch    = SAMTOOLS_COMPRESS_AND_SORT(aligned_bam_ch)
        indexed_bam_ch              = SAMTOOLS_INDEXING(compressed_and_sorted_ch)        

        _vcf_ch = GATK_HAPLOTYPECALLER(
            indexed_bam_ch,
            reference_fasta_ch,
            reference_fai_ch,
            reference_dict_ch
        )
    
    emit:
        _vcf_ch
        indexed_bam_ch
}