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