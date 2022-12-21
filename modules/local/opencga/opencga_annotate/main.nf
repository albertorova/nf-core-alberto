process OPENCGA_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'
    
    container 'tgbytes/opencga-base:2.4.5'
    containerOptions '--network host'

    input:
    tuple val(meta), path(vcf)

    //Plantear alguna manera de introducir los datos de forma dinamica

    shell:
    def user="arodriguez"
    def pass="Tsystems2022_arodriguez"
    def project="arodriguez@id_project_arodriguez"
    def study="arodriguez@id_project_arodriguez:id_study_arodriguez"

    """
    echo Acceder y linkar los archivos a OPENCGA

    /opt/opencga/bin/opencga.sh users login -u $user -p $pass

    /opt/opencga/bin/opencga.sh files create --study $study --path 'vcf' --type 'DIRECTORY'

    /opt/opencga/bin/opencga.sh files link -i /var/lib/mongo/tmp/alberto/results/variant_calling/haplotypecaller/joint_variant_calling/$vcf --path 'vcf' --study $study


    echo Indexado

    /opt/opencga/bin/opencga.sh operations variant-index --file $vcf --merge advanced
    

    echo Anotacion

    /opt/opencga/bin/opencga.sh operations variant-annotation-index --project $project

    /opt/opencga/bin/opencga.sh operations variant-stats-index --study $study --cohort 'ALL'

    /opt/opencga/bin/opencga.sh operations variant-secondary-index --project $project

    """
}
