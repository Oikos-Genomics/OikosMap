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
    Input individuals           : ${params.indlist}
    Directory of readfiles      : ${params.indir}
    Input refseq                : ${params.refseq}
    Number of threads           : ${params.threads}
    Output prefix               : ${params.prefix}
    """
    .stripIndent()

include { PRINT_HELP; CHECK_PARAMS_FOR_NULL; CHECK_FILE_FOR_EXISTENCE } from './modules/housekeeping_processes.nf'
include { CHECK_REFSEQ_FOR_INDEX; BWA_INDEX } from './modules/mapping_processes.nf'

help_message='''
    Basic Usage:
    nextflow run OikosMap.nf --indlist <names_of_inds.txt> --indir </path/to/directory/with/reads/> --refseq <refseq.fa> --prefix <output_prefix>
    
    Options:
        --help          Flag. Show this help message and exit
        --indlist       Mandatory
        --indir         Mandatory
        --refseq        Mandatory
        --threads       Defaults to 1/2 number on host machine
        --prefix        Name of output; defaults to 'out'. Should not contain whitespace.

    Notes:
        Nearly all options are mandatory.
        Any errors (of which there are probably many) should be reported to https://github.com/BirdmanRidesAgain/OikosMap/issues.
    '''

workflow {
    // INTRODUCTORY BEHAVIOR
    if (params.help) {
        PRINT_HELP(help_message)
        exit 10
    }

    // PARSE INPUTS
    CHECK_PARAMS_FOR_NULL([params.indlist, params.indir, params.refseq])
    CHECK_FILE_FOR_EXISTENCE([params.indlist, params.indir, params.refseq])
    CHECK_REFSEQ_FOR_INDEX(params.refseq) //FIXME - throws an error if the refseq isn't indexed. But we already index again inside NF, so redundant.

    //Add a backslash to --indir if there's not one
    clean_indir = ""
    if ( params.indir.substring(params.indir.length() - 1, params.indir.length()) != '/' ) { 
        clean_indir=params.indir+'/'
        } else { clean_indir = params.indir }

    fq_patterns = [
        clean_indir+'*{R1,R2}.fq.gz',
        clean_indir+'*{1,2}.fq.gz',
        clean_indir+'*{R1,R2}.fastq.gz',
        clean_indir+'*{1,2}.fastq.gz',
        clean_indir+'*{R1,R2}.fq',
        clean_indir+'*{1,2}.fq',
        clean_indir+'*{R1,R2}.fastq',
        clean_indir+'*{1,2}.fastq',
    ]

    reads_ch = Channel.fromFilePairs(fq_patterns)

    // CHECK FOR INDEXED REFSEQ
    refseq_ch = Channel.fromPath(params.refseq)
    BWA_INDEX(refseq_ch)
}
