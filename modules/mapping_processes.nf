/*
 * Keiler Collier
 * Helper script for OikosMap.nf, containing processes related to mapping.
 */


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

process BWA_INDEX {
    tag "Indexing ${refseq}"
    conda 'bwa'

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

process FASTP {
    tag "Trimming input ${sample_id}"
    publishDir "${params.prefix}_out/individuals/${sample_id}/fastp", mode: 'symlink', overwrite: 'false'
    conda 'fastp'
    maxForks params.threads

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
    publishDir "${params.prefix}_out/individuals/${sample_id}/mapping", mode: 'symlink', overwrite: 'false'
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

process VCF_CALL_MASS {
    tag "Variant-calling all mapped files"
    publishDir "${params.prefix}_out/vcf", mode: 'symlink', overwrite: 'false'
    conda 'bcftools'

    input:
    path(bam)
    tuple path(refseq), path(amb), path(ann), path(btw), path(pac), path(sa)
    val(prefix)

    output:
    tuple val(prefix), path("${prefix}.vcf.gz"), emit: 'vcf_mass'

    script:
    """
    FORMAT_STRING=FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR
    find . -name "*final.bam" | xargs bcftools mpileup -f ${refseq} -q 60 -a "\$FORMAT_STRING" -Ou | bcftools call -m -a GQ,GP --variants-only -Oz -o "${prefix}.vcf.gz"
    """
}

process VCF_CALL_IND {
    tag "Variant-calling ${sample_id}"
    publishDir "${params.prefix}_out/individuals/${sample_id}/mapping", mode: 'symlink', overwrite: 'false'
    conda 'bcftools'
    
    input:
    tuple val(sample_id), path(bam), path(refseq), path(amb), path(ann), path(btw), path(pac), path(sa)

    output:
    path("${sample_id}.vcf.gz"), emit: 'vcf_ind'

    script:
    """
    FORMAT_STRING=FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR
    SAMPLE_ID=\$(basename -s _final.bam ${bam})
    find . -name "${bam}" | xargs bcftools mpileup -f ${refseq} -q 60 -a "\$FORMAT_STRING" -Ou | bcftools call -m -a GQ,GP --variants-only -Oz -o "\${SAMPLE_ID}.vcf.gz"
    """
}

