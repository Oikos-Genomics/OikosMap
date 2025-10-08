/*
 * Keiler Collier
 * OikosMap V1 - used to map and variant-call shortreads.
 *
 * Started: 12 Sep 2025
 * Last update: 12 Sep 2025
 * Pipeline input params supplied from nextflow.config
 */

log.info """\
    O I K O S M A P - N F   P I P E L I N E
    ===================================
    Project directory           : $projectDir
    Directory of readfiles      : ${params.indir}
    Input refseq                : ${params.refseq}
    Number of threads           : ${params.threads}
    Output prefix               : ${params.prefix}
    """
    .stripIndent()

include { PRINT_HELP; CHECK_PARAMS_FOR_NULL; CHECK_FILE_FOR_EXISTENCE } from './modules/housekeeping_processes.nf'
include { CHECK_REFSEQ_FOR_INDEX; BWA_INDEX; FASTP; BWA_MEM; VCF_CALL_MASS; VCF_CALL_IND } from './modules/mapping_processes.nf'

help_message='''
    Basic Usage:
    nextflow run OikosMap.nf --indir </path/to/directory/with/reads/> --suffix <_{R1,R2}.fq.gz> --refseq <refseq.fa> --prefix <output_prefix>
    
    Options:
        --help          Optional. Show this help message and exit
        --indir         Mandatory. The directory where all the read files are located
        --suffix        Mandatory. The suffix of all read files in --indir. Defaults to '_{R1,R2}.fq.gz'
        --refseq        Mandatory. Your reference sequence
        --ind_vcfs      Optional. Set to generate single vcfs per individual. May increase speed.
        --threads       Optional. Defaults to 1/2 number on host machine
        --prefix        Optional. The name of the output directory. Defaults to 'out'.

    Notes:
        Any errors (of which there are many) should be reported to https://github.com/BirdmanRidesAgain/OikosMap/issues.
    '''

//nextflow.preview.output = true

workflow {
    //Main workflow verifies inputs to ensure quality
    main:
    // INTRODUCTORY BEHAVIOR
    if (params.help) {
        PRINT_HELP(help_message)
        exit 10
    }

    // PARSE INPUTS
    CHECK_PARAMS_FOR_NULL([params.indir, params.refseq])
    CHECK_FILE_FOR_EXISTENCE([params.indir, params.refseq])
    //CHECK_REFSEQ_FOR_INDEX(params.refseq) //FIXME - throws an error if the refseq isn't indexed. But we already index again inside NF, so redundant.
    //Add a backslash to --indir if there's not one
    clean_indir = ""
    if ( params.indir.substring(params.indir.length() - 1, params.indir.length()) != '/' ) { 
        clean_indir=params.indir+'/'
        } else { clean_indir = params.indir }

    fq_patterns = [
        clean_indir+'**'+params.suffix,
    ]

    reads_ch = Channel.fromFilePairs(fq_patterns, flat: true)
    refseq_ch = Channel.fromPath(params.refseq)

    //MAP_AND_VARCALL(reads_ch, refseq_ch)
    BWA_INDEX(refseq_ch)
    FASTP(reads_ch)

    BWA_MEM(FASTP.out.fq_trimmed.combine(BWA_INDEX.out.indexed_refseq), params.threads)

    if ( params.ind_vcfs ) {
        VCF_CALL_IND(BWA_MEM.out.bamfile_ind.combine(BWA_INDEX.out.indexed_refseq))
    } else { VCF_CALL_MASS(BWA_MEM.out.bamfile_mass.collect(), BWA_INDEX.out.indexed_refseq, params.prefix)}
    
    //publish:
    //fastp_QC = fastp

}



//BWA_MEM.out.bam_ch.collect()
/*
output {
    fastp_QC {

    }
}

workflow MAP_AND_VARCALL {
    take:
    reads
    refseq

    main:
    BWA_INDEX(refseq)
    FASTP(reads)

    emit:
    fastp = FASTP.out.fq_QC
}
*/
