process FREEBAYES {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::freebayes=1.3.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.6--hbfe0e7f_2' :
        'quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2' }"

    input:
    tuple val(meta), path(input_1)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input            = input_1        ? "${input_1}"        : "${input_1}"

    """
    freebayes \\
        -f $fasta \\
        $args \\
        $input > ${prefix}.vcf

    bgzip ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}
