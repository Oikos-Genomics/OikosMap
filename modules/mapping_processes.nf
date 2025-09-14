/*
 * Keiler Collier
 * Helper script for OikosMap.nf, containing processes related to mapping.
 */


def CHECK_REFSEQ_FOR_INDEX(refseq) {
    // CHECK FOR INDEXED REFSEQ
    refseq_ch = Channel.fromPath([refseq,refseq+'.{amb,ann,bwt,pac,sa}'])
    refseq_ch.count().map { n ->
    if( n < 6 ) {
        error("${refseq} is not indexed. Index with `bwa mem ${refseq}`, and rerun.")
        }
    }
}

process BWA_INDEX {
    tag "Indexing ${refseq}"    
    conda 'bwa'

    input:
    path refseq
    
    output:
    stdout emit: bwa_check
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


//process FASTP {
//    tag "Trimming input reads"
//    cpus 10
//    conda 'bioconda::fastp'
//
//    //input:
//    //tuple val(sample_id), val(trim_status), path(reads)
//    //
//    //output:
//    //path "${sample_id}_${trim_status}_NanoPlot"
//
//    script:
//    println("hello world")
//}
//
//process BWA_MEM {
//    tag "Trimming input reads"
//    cpus 10
//    conda 'bioconda::fastp'
//
//    //input:
//    //tuple val(sample_id), val(trim_status), path(reads)
//    //
//    //output:
//    //path "${sample_id}_${trim_status}_NanoPlot"
//
//    script:
//    println("hello world")
//}
//
//process BCFTOOLS {
//    tag "Trimming input reads"
//    cpus 10
//    conda 'bioconda::fastp'
//
//    //input:
//    //tuple val(sample_id), val(trim_status), path(reads)
//    //
//    //output:
//    //path "${sample_id}_${trim_status}_NanoPlot"
//
//    script:
//    println("hello world")
//}

