process DO_BWA_MEM {
    tag "Mapping input ${sample_id}"
    publishDir "${params.prefix}_out/individuals/${sample_id}/mapping", mode: 'copy', overwrite: 'false'
    conda 'bwa samtools'
    maxForks 1

    input:
    tuple val(sample_id), path(read_1), path(read_2), path(refseq), path(amb), path(ann), path(btw), path(pac), path(sa)
    val(threads)

    output:
    path("${sample_id}_final.bam"), emit: bamfile_mass
    tuple val(sample_id), path("${sample_id}_final.bam"), emit: bamfile_ind
    

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