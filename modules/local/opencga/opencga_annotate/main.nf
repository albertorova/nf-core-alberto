process OPENCGA_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'
    
    //Imagen Docker de Opencga
    container 'tgbytes/opencga-base:2.4.5'
    containerOptions '--network host'

    input:
    tuple val(meta), path(vcf)

    //Plantear en el futuro una forma de incluir datos dinamicamente

    shell:
    def user="arodriguez"
    def pass="Tsystems2022_arodriguez"
    def project="arodriguez@id_project_arodriguez"
    def study="arodriguez@id_project_arodriguez:id_study_arodriguez"

    """
    echo Acceder y linkar los archivos a OPENCGA

    opencga.sh users login -u $user -p $pass

    opencga.sh files link -i /home/albertorova/nextflow/pruebas/results/strelka/$vcf --study $study


    echo Indexado

    opencga.sh operations variant-index --file $vcf
    

    echo Anotacion

    opencga.sh operations variant-annotation-index --project $project

    opencga.sh operations variant-stats-index --study $study --cohort 'ALL'

    opencga.sh operations variant-secondary-index --project $project

    """
}
