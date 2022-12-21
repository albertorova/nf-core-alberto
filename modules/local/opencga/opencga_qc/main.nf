process OPENCGA_QC {
    tag "$meta.id"
    label 'process_medium'
    
    container 'tgbytes/opencga-base:2.4.5'
    containerOptions '--network host'

    input:
    tuple val(meta), path(bam)
    path(bai)
    tuple val(meta2), path(bed)
    path(dict)

    //Plantear alguna manera de introducir los datos de forma dinamica

    shell:
    def user="arodriguez"
    def pass="Tsystems2022_arodriguez"
    def project="arodriguez@id_project_arodriguez"
    def study="arodriguez@id_project_arodriguez:id_study_arodriguez"
    def individual="ISDBM322015"
    def sample="ISDBM322015"
    def family="corpas"
    def genes="DDX11L1"


    """
    echo $bam
    echo $bai
    echo $bed
    echo $dict

    echo Acceder y linkar los archivos a OPENCGA

    /opt/opencga/bin/opencga.sh users login -u $user -p $pass

    /opt/opencga/bin/opencga.sh files create --study $study --path 'data' --type 'DIRECTORY'

    /opt/opencga/bin/opencga.sh files link -i /var/lib/mongo/nextflow/pruebas/results/preprocessing/recalibrated/$meta.id/$bam --study $study --path 'data'

    /opt/opencga/bin/opencga.sh files link -i /var/lib/mongo/nextflow/pruebas/results/preprocessing/recalibrated/$meta.id/$bai --study $study --path 'data'

    /opt/opencga/bin/opencga.sh files link -i /var/lib/mongo/nextflow/pruebas/results/bedtools/$bed --study $study --path 'data'

    /opt/opencga/bin/opencga.sh files link -i /var/lib/mongo/nextflow/pruebas/results/samtools/$dict --study $study --path 'data'

    
    echo QC VARIANT STATS

    /opt/opencga/bin/opencga.sh operations variant-stats-index --study $study --cohort 'ALL'


    echo QC Sex Inference/Mendelian Errors

    /opt/opencga/bin/opencga.sh variant individual-qc-run --individual $individual --sample $sample --study $study


    echo QC Relatedness

    /opt/opencga/bin/opencga.sh variant family-qc-run --family $family --study $study


    echo QC Samtools Stats/Plots/Flagstats
    
    /opt/opencga/bin/opencga.sh alignments qc-run --bam-file $bam --bed-file $bed --dict-file $dict --study $study


    echo QC Gene Coverage Stats

    /opt/opencga/bin/opencga.sh alignments coverage-qc-genecoveragestats-run --bam-file $bam --genes $genes --study $study


    echo QC Sample stats

    /opt/opencga/bin/opencga.sh variant sample-qc-run --sample $sample --study $study

    """

}
