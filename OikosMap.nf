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

def print_help() {
    //Prints out OikosMap instructions.
    log.info"""
    Basic Usage:
    nextflow run OikosMap.nf [--ONT_raw OR --PB_raw]
    
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
    
    """.stripIndent()
}

include { INDEX_REFSEQ } from './modules/mapping_processes.nf'

workflow {
    // INTRODUCTORY BEHAVIOR
    if (params.help) {
        print_help()
        exit 0
    }

    //CHECK INPUT PARAMETERS:
    if ( params.indlist == null | params.indir == null | params.refseq == null ) {
        def input_errors = []
        if (params.indlist == null) {
            input_errors.add("\n\t--indlist not provided")
        }

        if (params.indir == null) {
            input_errors.add("--indir not provided")
        }

        if (params.refseq == null) {
            input_errors.add("--refseq not provided")
        }

        if (input_errors.size() > 0) {
            error("Input files missing: ${input_errors.join('\n\t')}\nSee --help for usage instructions.")
        }
    }

    // PARSE INPUTS

    //check for backslash and existence on --indir
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
    reads_ch.view()

    INDEX_REFSEQ()

}
