/*
 * Keiler Collier
 * OikosMap V1 - used to map and variant-call shortreads.
 *
 * Started: 12 Sep 2025
 * Last update: 12 Sep 2025
 * Pipeline input params supplied from nextflow.config
 */

    help = false
    indlist = null
    indir = null
    refseq = null 
    threads = max_cpus/2
    prefix = out


log.info """\
    O I K O S M A P - N F   P I P E L I N E
    ===================================
    Project directory           : $projectDir
    Input individuals           : ${params.indlist}
    Directory of readfiles      : ${params.indir}
    Input refseq                : ${params.refseq}
    Number of threads           : ${params.threads}
    Output prefix               : ${params.out}
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

//include { REMOVE_CONTAMINANTS } from './modules/processes.nf' //initial filtering for contaminants
//include { NANOPLOT as NANOPLOT_RAW; NANOPLOT as NANOPLOT_TRIM; NANOFILT;  //nanoplot-related
//} from './modules/mapping_processes.nf'
//include { VARIANT_MPILEUP_CALL } from './modules/variant_processes.nf' //merged assembly-related

workflow {
// INTRODUCTORY BEHAVIOR
    if (params.help) {
        print_help()
        exit 0
    }
    //CHECK INPUT PARAMETERS:
    if ( params.indlist == null | params.indir == null | params.refseq == null ) {
        error "One or more missing arguments. See --help for usage instuctions." 
    }
}
