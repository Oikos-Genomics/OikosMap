#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Keiler Collier
 * OikosMap V1 - used to map and variant-call shortreads.
 *
 * Started: 12 Sep 2025
 * Last update: 12 Sep 2025
 * Pipeline input params supplied from nextflow.config
 */

// import modules
include { PRINT_HELP } from './modules/print_help'
include { CHECK_PARAMS_FOR_NULL } from './modules/check_params_for_null'
include { CHECK_FILE_FOR_EXISTENCE } from './modules/check_file_for_existence'
include { CHECK_REFSEQ_FOR_INDEX } from './modules/check_refseq_for_index'
include { DO_BWA_INDEX } from './modules/do_bwa_index'
include { DO_FASTP } from './modules/do_fastp'
include { DO_BWA_MEM } from './modules/do_bwa_mem'
include { DO_VCF_CALL_MASS } from './modules/do_vcf_call_mass'
include { DO_VCF_CALL_IND } from './modules/do_vcf_call_ind'


//nextflow.preview.output = true

workflow {
    //Main workflow verifies inputs to ensure quality
    main:
    // INTRODUCTORY BEHAVIOR
    if (params.help) {
        PRINT_HELP()
        exit 10
    }

    log.info """\
    O I K O S M A P - N F   P I P E L I N E
    ===================================
    Project directory           : $projectDir
    Directory of readfiles      : ${params.indir}
    Readfile suffix             : ${params.suffix}
    Input refseq                : ${params.refseq}
    Individual VCFs             : ${params.ind_vcfs}
    Number of threads           : ${params.threads}
    Output prefix               : ${params.prefix}
    """
    .stripIndent()

    // PARSE INPUTS
    CHECK_PARAMS_FOR_NULL([params.indir, params.refseq])
    CHECK_FILE_FOR_EXISTENCE([params.indir, params.refseq])
    IS_REFSEQ_INDEXED=CHECK_REFSEQ_FOR_INDEX(params.refseq)

    if ( IS_REFSEQ_INDEXED ) {
        indexed_refseq_ch=channel.fromPath([params.refseq,params.refseq+'.{amb,ann,bwt,pac,sa}']).collect()
        } else {
            indexed_refseq_ch=DO_BWA_INDEX(channel.fromPath(params.refseq))
        }

    //Add a backslash to --indir if there's not one
    clean_indir = ""
    if ( params.indir.substring(params.indir.length() - 1, params.indir.length()) != '/' ) { 
        clean_indir=params.indir+'/'
        } else { clean_indir = params.indir }
    fq_patterns = [
        clean_indir+'**'+params.suffix,
    ]
    reads_ch = channel.fromFilePairs(fq_patterns, flat: true)

    // TRIM, MAP, VARIANT-CALL
    DO_FASTP(reads_ch)
    DO_BWA_MEM(DO_FASTP.out.fq_trimmed.combine(indexed_refseq_ch), params.threads)

    if ( params.ind_vcfs ) {
        DO_VCF_CALL_IND(DO_BWA_MEM.out.bamfile_ind.combine(indexed_refseq_ch))
    } else { DO_VCF_CALL_MASS(DO_BWA_MEM.out.bamfile_mass.collect(), indexed_refseq_ch, params.prefix)}

   //publish:
   //fastp_QC = fastp

/*
workflow MAP_AND_VARCALL {
    take:
    reads
    refseq

    main:
    BWA_INDEX(refseq)
    FASTP(reads)

    emit:
    fastp = FASTP.out.fq_QC
*/
}

