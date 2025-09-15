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
    //stdout emit: bwa_check
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

process FASTP {
    tag "Trimming input ${sample_id}"
    publishDir "${params.prefix}_out/${sample_id}/fastp", mode: 'symlink', overwrite: 'false'
    conda 'fastp'

    input:
    tuple val(sample_id), path(read_1), path(read_2)
    
    output:
    tuple val(sample_id), path("${read_1_out}.fq.gz"), path("${read_2_out}.fq.gz"), emit: fq_trimmed
    tuple path("${sample_id}.html"), path("${sample_id}.json"), emit: fq_QC

    script:
    read_1_out = sample_id+"_fastp_1"
    read_2_out = sample_id+"_fastp_2"
    """
    fastp -i ${read_1} -I ${read_2} -o ${read_1_out}.fq.gz -O ${read_2_out}.fq.gz
    mv fastp.html ${sample_id}.html
    mv fastp.json ${sample_id}.json
    """
}

process BWA_MEM {
    tag "Mapping input ${sample_id}"
    publishDir "${params.prefix}_out/${sample_id}/mapping", mode: 'symlink', overwrite: 'false'
    conda 'bwa samtools'

    input:
    tuple val(sample_id), path(read_1), path(read_2), path(refseq), path(amb), path(ann), path(btw), path(pac), path(sa), val(threads)
    val(threads)

    output:
    tuple val(sample_id), path("${sample_id}_final.bam"), emit: bamfile

    script:
    """
    bwa mem -t "${threads}" "${refseq}" "${read_1}" "${read_2}" | samtools view -@ "${threads}" -bS > "${sample_id}.bam"

    samtools collate -@ "${threads}" -o ${sample_id}_namecollate.bam ${sample_id}.bam
    samtools fixmate -@ "${threads}" -m ${sample_id}_namecollate.bam ${sample_id}_fixmate.bam
    samtools sort -@ "${threads}" -o ${sample_id}_positionsort.bam ${sample_id}_fixmate.bam
    samtools markdup -r -@ "${threads}" ${sample_id}_positionsort.bam ${sample_id}_markdup.bam
    samtools addreplacerg -@ "${threads}" -r ID:${sample_id} -r LB:${sample_id} -r SM:${sample_id} -o ${sample_id}_final.bam ${sample_id}_markdup.bam

    rm ${sample_id}_{namecollate,fixmate,positionsort,markdup}.bam 
    rm ${sample_id}.bam
    """
}
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

