process DO_BWA_INDEX {
    tag "Indexing ${refseq}"

    input:
    path refseq
    
    output:
    tuple path(refseq), path("*amb"), path("*ann"), path("*bwt"), path("*pac"), path("*sa"), emit: indexed_refseq

    script:

    """ 
    # Check if BWA index exists
    if [ ! -f "${refseq}.amb" ]; then
        echo "Creating BWA index..."
        bwa index ${refseq}
    else
        echo "Found existing BWA index files"
    fi
    """
}