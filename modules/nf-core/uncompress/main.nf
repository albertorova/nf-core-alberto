process UNCOMPRESS {

    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(fastq_gz)

    output:
    //path fastq_gz.baseName
    tuple val(meta), path("*.vcf"), emit: uncompress_vcf

    script:
    """
    gunzip -c "${fastq_gz}" > "${fastq_gz.baseName}"
    """
}