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