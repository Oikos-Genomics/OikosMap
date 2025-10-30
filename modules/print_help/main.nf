def PRINT_HELP() {
def message='''
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
    log.info message.stripIndent()
}