process DO_FASTP {
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