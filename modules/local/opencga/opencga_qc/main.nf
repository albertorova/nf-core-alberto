process OPENCGA_QC {
    tag "$meta.id"
    label 'process_medium'
    
    //Imagen Docker de Opencga
    container 'tgbytes/opencga-base:2.4.5'
    containerOptions '--network host'

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bed)
    path(dict)

    //Plantear en el futuro una forma de incluir datos dinamicamente

    shell:
    def user="arodriguez"
    def pass="Tsystems2022_arodriguez"
    def project="arodriguez@id_project_arodriguez"
    def study="arodriguez@id_project_arodriguez:id_study_arodriguez"
    def family="corpas"


    """
    echo $bam
    echo $bed
    echo $dict

    echo Acceder y linkar los archivos a OPENCGA

    opencga.sh users login -u $user -p $pass

    opencga.sh files link -i /home/albertorova/nextflow/pruebas/results/gatk4/$bam --study $study

    opencga.sh files link -i /home/albertorova/nextflow/pruebas/results/bedtools/$bed --study $study

    opencga.sh files link -i /home/albertorova/nextflow/pruebas/results/gatk4/$dict --study $study

    
    echo QC VARIANT STATS

    opencga.sh operations variant-stats-index --study $study --cohort 'ALL'


    echo QC Relatedness

    opencga.sh variant family-qc-run --family $family --study $study


    echo QC Samtools Stats/Plots/Flagstats
    
    opencga.sh alignments qc-run --bam-file $bam --bed-file $bed --dict-file $dict --study $study

    """

}
