process UNCOMPRESS {

    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(archive)

   // output:
   // tuple val(meta), path("${archive.baseName}/"), emit: unzipped_archive

    script:
    def a = PWD
    def b = "/results/freebayes/"
    """
    gunzip $a$b$archive

    """
}