nextflow.enable.dsl=2

include { FASTQC }          from '../../modules/fastqc/main.nf'
include { MULTIQC }         from '../../modules/multiqc/main.nf'

workflow QUALITY_CONTROL {
    take:
        reads_ch

    main:
        fastqc_reports_ch = FASTQC(reads_ch)
        MULTIQC(fastqc_reports_ch.collect())
    
    emit:
        fastqc_reports_ch

}