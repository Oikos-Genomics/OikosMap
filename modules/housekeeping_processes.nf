/*
 * Keiler Collier
 * Helper script for OikosMap.nf, containing processes related to program inputs/outputs/options.
 */

def PRINT_HELP(message) {
    log.info message.stripIndent()
}

def CHECK_PARAMS_FOR_NULL (inlist) {
    //checks that each item is non-null
    inlist.each { i -> 
        def in_errors = ["\n"]
        if ( "${i}" == null ) {
            in_errors.add("${i} not provided")
        }
        if ( in_errors.size() > 1 ) {
            error("Input params missing: ${in_errors.join('\n\t')}\nSee --help for usage instructions.")
        }
    }
}

def CHECK_FILE_FOR_EXISTENCE (inlist) {
    //checks that each file exists
    inlist.each { i -> 
        def in_errors = ["\n"]
        if (!file(i).exists()) {
            in_errors.add("${i} does not exist")
        }
        if ( in_errors.size() > 1 ) {
            error("Input files missing: ${in_errors.join('\n\t')}\nSee --help for usage instructions.")
        }
    }
}

