/*
 * Keiler Collier
 * Helper script for OikosMap.nf, containing processes related to mapping.
 */

process INDEX_REFSEQ {
    tag "Indexing refseq with bwa index"
    cpus 10
    conda 'bioconda::bwa'

    //input:
    //tuple val(sample_id), val(trim_status), path(reads)
    //
    //output:
    //path "${sample_id}_${trim_status}_NanoPlot"

    script:
    println("hello world")

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

