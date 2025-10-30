process DO_VCF_CALL_IND {
    tag "Variant-calling ${sample_id}"
    
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