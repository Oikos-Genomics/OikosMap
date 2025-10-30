def CHECK_REFSEQ_FOR_INDEX(refseq) {
    def indexFiles = [
        new File(refseq),
        new File(refseq + '.amb'),
        new File(refseq + '.ann'), 
        new File(refseq + '.bwt'),
        new File(refseq + '.pac'),
        new File(refseq + '.sa')
    ]
    
    return indexFiles.every { it.exists() }
}
