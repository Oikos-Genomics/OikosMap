process DO_VCF_CALL_MASS {
    tag "Variant-calling all mapped files"
    publishDir "${params.prefix}_out/vcf", mode: 'move', overwrite: 'false'
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