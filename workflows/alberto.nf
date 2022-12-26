/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowAlberto.initialise(params, log)


// Check de los parametros
def checkPathParamList = [ params.input, 
                           params.multiqc_config, 
                           params.dbsnp,
                           params.fasta,
                           params.intervals,
                           params.known_indels,
                           params.snpeff_cache,
                           params.vep_cache,
                           params.germline_resource,
                            ]



for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules

include { INPUT_CHECK } from '../subworkflows/local/input_check'

dbsnp              = params.dbsnp                  ? Channel.fromPath(params.dbsnp).collect()                        : Channel.empty()
fasta              = params.fasta                  ? Channel.fromPath(params.fasta).collect()                        : Channel.empty()
intervals          = params.intervals              ? Channel.fromPath(params.intervals).collect()                    : Channel.empty()
known_indels       = params.known_indels           ? Channel.fromPath(params.known_indels).collect()                 : Channel.empty()

snpeff_db          = params.snpeff_db              ?: Channel.empty()
vep_cache_version  = params.vep_cache_version      ?: Channel.empty()

snpeff_cache       = params.snpeff_cache           ? Channel.fromPath(params.snpeff_cache).collect()                 : []
vep_cache          = params.vep_cache              ? Channel.fromPath(params.vep_cache).collect()                    : []

germline_resource  = params.germline_resource      ? Channel.fromPath(params.germline_resource).collect()            : Channel.value([]) 


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Modulos nf-core/modules


include { GATK4_CREATESEQUENCEDICTIONARY         } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                         } from '../modules/nf-core/samtools/faidx/main'

include { BWA_INDEX as BWAMEM1_INDEX             } from '../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX                          } from '../modules/nf-core/bwamem2/index/main'

include { BOWTIE2_BUILD                          } from '../modules/nf-core/bowtie2/build/main'

include { TABIX_TABIX as TABIX_DBSNP             } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS      } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE } from '../modules/nf-core/tabix/tabix/main'

include { FASTQC                                 } from '../modules/nf-core/fastqc/main'

include { TRIMMOMATIC                            } from '../modules/nf-core/trimmomatic/main'
include { FASTP                                  } from '../modules/nf-core/fastp/main'
include { TRIMGALORE                             } from '../modules/nf-core/trimgalore/main'

include { BWA_MEM as BWAMEM1_MEM                 } from '../modules/nf-core/bwa/mem/main'
include { BWAMEM2_MEM                            } from '../modules/nf-core/bwamem2/mem/main'

include { BOWTIE2_ALIGN                          } from '../modules/nf-core/bowtie2/align/main'

include { SAMBLASTER                             } from '../modules/nf-core/samblaster/main'

include { FGBIO_GROUPREADSBYUMI                  } from '../modules/nf-core/fgbio/groupreadsbyumi/main'

include { GATK4_MARKDUPLICATES_SPARK             } from '../modules/nf-core/gatk4/markduplicatesspark/main'
include { GATK4_ESTIMATELIBRARYCOMPLEXITY        } from '../modules/nf-core/gatk4/estimatelibrarycomplexity/main'

include { PICARD_ADDORREPLACEREADGROUPS          } from '../modules/nf-core/picard/addorreplacereadgroups/main'

include { GATK4_BASERECALIBRATOR                 } from '../modules/nf-core/gatk4/baserecalibrator/main'

include { GATK4_APPLYBQSR                        } from '../modules/nf-core/gatk4/applybqsr/main'

include { SAMTOOLS_INDEX                         } from '../modules/nf-core/samtools/index/main'

include { BEDTOOLS_BAMTOBED                      } from '../modules/nf-core/bedtools/bamtobed/main'

include { FREEBAYES                              } from '../modules/nf-core/freebayes/main'
include { STRELKA_GERMLINE                       } from '../modules/nf-core/strelka/germline/main'
include { GATK4_MUTECT2                          } from '../modules/nf-core/gatk4/mutect2/main'

include { SNPEFF                                 } from '../modules/nf-core/snpeff/main'

include { UNCOMPRESS                             } from '../modules/nf-core/uncompress/main'

include { VCF2MAF                                } from '../modules/nf-core/vcf2maf/main'

include { BCFTOOLS_STATS                         } from '../modules/nf-core/bcftools/stats/main' 
include { VCFTOOLS                               } from '../modules/nf-core/vcftools/main'
include { RTGTOOLS_VCFEVAL                       } from '../modules/nf-core/rtgtools/vcfeval/main'

include { OPENCGA_QC                             } from '../modules/local/opencga/opencga_qc/main'

include { OPENCGA_ANNOTATE                       } from '../modules/local/opencga/opencga_annotate/main'

include { MULTIQC                                } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS            } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ALBERTO {

    //Channels para reports de diferentes herramientas

        ch_reports  = Channel.empty()

        ch_versions = Channel.empty()

        qc_reports  = Channel.empty()

        ch_logs = Channel.empty()

//*****************************************************************************************************************

    //Preparar fasta, fai, dict
    
        GATK4_CREATESEQUENCEDICTIONARY(fasta)
        SAMTOOLS_FAIDX(fasta.map{ it -> [[id:it[0].baseName], it] }) 
    
        dict                             = GATK4_CREATESEQUENCEDICTIONARY.out.dict                               // path: genome.fasta.dict
        fasta_fai                        = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }                      // path: genome.fasta.fai
        
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

//*****************************************************************************************************************

    //if (params.mapeador == 'bwamem') {

    //Preparar indices para BWAMEM/BWAMEM2

        //BWAMEM1_INDEX(fasta.map{ it -> [[id:it[0].baseName], it] }) // If aligner is bwa-mem
        //BWAMEM2_INDEX(fasta.map{ it -> [[id:it[0].baseName], it] }) // If aligner is bwa-mem2

        //bwa                              = BWAMEM1_INDEX.out.index.map{ meta, index -> [index] }.collect()       // path: bwa/*
        //bwamem2                          = BWAMEM2_INDEX.out.index.map{ meta, index -> [index] }.collect()       // path: bwamem2/*

        //ch_versions = ch_versions.mix(BWAMEM1_INDEX.out.versions)
        //ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

//*****************************************************************************************************************

    //} else {

    //Preparar indices para BOWTIE2

        BOWTIE2_BUILD(fasta.map{ it -> [[id:it[0].baseName], it] }) // If aligner is bowtie2

        bowtie2                          = BOWTIE2_BUILD.out.index.map{ meta, index -> [index] }.collect()       // path: bwamem2/*

        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)

    //}
    
//*****************************************************************************************************************

    //Preparar KNOWN_INDELS y DBSNP

        TABIX_KNOWN_INDELS( known_indels.flatten().map{ it -> [[id:it.baseName], it] } )
        TABIX_DBSNP(dbsnp.flatten().map{ it -> [[id:it.baseName], it] })
        TABIX_GERMLINE_RESOURCE(germline_resource.flatten().map{ it -> [[id:it.baseName], it] })

        known_indels_tbi                 = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [tbi] }.collect()        // path: {known_indels*}.vcf.gz.tbi
        dbsnp_tbi                        = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }.collect()               // path: dbsnb.vcf.gz.tbi
        germline_resource_tbi            = TABIX_GERMLINE_RESOURCE.out.tbi.map{ meta, tbi -> [tbi] }.collect()   // path: germline_resource.vcf.gz.tbi

        ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)
        ch_versions = ch_versions.mix(TABIX_KNOWN_INDELS.out.versions)
        ch_versions = ch_versions.mix(TABIX_GERMLINE_RESOURCE.out.versions)

        known_sites_indels     = dbsnp.concat(known_indels).collect()
        known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()



//*****************************************************************************************************************
    
    //
    //VALIDAR INPUT
    //
        INPUT_CHECK(ch_input)
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

//*****************************************************************************************************************

    //
    //FASTQC
    //
        //FASTQC (INPUT_CHECK.out.reads)

        //ch_reports  = ch_reports.mix(FASTQC.out.zip.collect{meta, logs -> logs})
        //ch_versions = ch_versions.mix(FASTQC.out.versions.first())

//*****************************************************************************************************************

    //FASTP/Trimming

        trimmed_reads  = Channel.empty()

        //save_trimmed_fail = false
        //save_merged = false
        //FASTP(INPUT_CHECK.out.reads,[],save_trimmed_fail,save_merged)

        //trimmed_reads  = FASTP.out.reads

        //ch_reports = ch_reports.mix(FASTP.out.json.collect{meta, json -> json},FASTP.out.html.collect{meta, html -> html})

//*****************************************************************************************************************

    //Trimmomatic SOLUCIONADO PROBLEMA:
    //AÑADIDA LINEA DE ARGUMENTOS --> ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

        TRIMMOMATIC(INPUT_CHECK.out.reads)

        trimmed_reads = trimmed_reads.mix(TRIMMOMATIC.out.trimmed_reads)

        ch_logs = ch_logs.mix(TRIMMOMATIC.out.log.first())
        ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())


//*****************************************************************************************************************

    //Trimgalore 
       
        //TRIMGALORE(INPUT_CHECK.out.reads) 

        //trimmed_reads = trimmed_reads.mix(TRIMGALORE.out.reads)

        //ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())


//*****************************************************************************************************************

    //if (params.mapeador == 'bwamem') {

    // Mapeado con BWAMEM/BWAMEM2    

        //sort_bam = true        
        //BWAMEM1_MEM(trimmed_reads,   bwa.map{ it -> [[id:it[0].baseName], it] }, sort_bam) // If aligner is bwa-mem
        //BWAMEM2_MEM(trimmed_reads,   bwamem2.map{ it -> [[id:it[0].baseName], it] }, sort_bam) // If aligner is bwa-mem2

        //Solo nos quedamos con el BWAMEM2
        //ch_bam_mapped = BWAMEM1_MEM.out.bam
        //ch_bam_mapped = BWAMEM2_MEM.out.bam

        //ch_versions = ch_versions.mix(BWAMEM1_MEM.out.versions.first())
        //ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())


    //} else {

    // Mapeado con BOWTIE2  

        save_unaligned = false     
        sort_bam = true   
        BOWTIE2_ALIGN(trimmed_reads,   bowtie2.map{ it -> [[id:it[0].baseName], it] },save_unaligned,sort_bam) 

        ch_bam_mapped = BOWTIE2_ALIGN.out.bam

        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    //}


//*****************************************************************************************************************
    
    //AÑADIR READ GROUPS AL BAM
    
        //SAMBLASTER(ch_bam_mapped)

        //FGBIO_GROUPREADSBYUMI
        //group_by_umi_strategy = 'Adjacency'
        //FGBIO_GROUPREADSBYUMI(SAMBLASTER.out.bam, group_by_umi_strategy)

        //PICARD_ADDORREPLACEREADGROUPS ARREGLAR!!
        //Añade grupos al bam. Tambien puede usarse para indexarlo
        //PICARD_ADDORREPLACEREADGROUPS(ch_bam_mapped)


//*****************************************************************************************************************

    //MARKDUPLICATES ARREGLAR!!

        //GATK4_MARKDUPLICATES_SPARK(ch_bam_mapped, fasta, fasta_fai, dict)

        //GATK4_ESTIMATELIBRARYCOMPLEXITY(ch_bam_mapped, fasta, fasta_fai, dict)

        //Reports
        //qc_reports = qc_reports.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.metrics)

    
        //ch_versions = ch_versions.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.versions.first())
        //ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES_SPARK.out.versions)


//******************************************************************************************************************
   
    //Indexar BAM

        SAMTOOLS_INDEX(ch_bam_mapped)

        ch_bam_bai = SAMTOOLS_INDEX.out.bai_solo

        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())


//*****************************************************************************************************************

    //OBTENER BED DEL BAM

        //BEDTOOLS_BAMTOBED(ch_bam_mapped)

        //ch_bed = BEDTOOLS_BAMTOBED.out.bed
        //ch_bed_solo = BEDTOOLS_BAMTOBED.out.bed_solo

        //ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions.first())

//*****************************************************************************************************************
       
    //GENERAR TABLAS PARA RECALIBRADO
        
        //PRIMERO GENERAMOS LAS TABLAS ARREGLAR!!
        //ERROR: A USER ERROR has occurred: Number of read groups must be >= 1, but is 0

        //GATK4_BASERECALIBRATOR(BWAMEM2_MEM.out.bam,ch_bam_bai,intervals,fasta,fasta_fai,dict,known_sites_indels,known_sites_indels_tbi)

        //ch_table = GATK4_BASERECALIBRATOR.out.table

        //ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)


//*****************************************************************************************************************

    //Recalibrado

        //GATK4_APPLYBQSR(BWAMEM2_MEM.out.bam,ch_bam_bai,ch_table,intervals,fasta,fasta_fai,dict)
        
        //recal_bam = GATK4_APPLYBQSR.out.bam

        //ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)


//*****************************************************************************************************************


    //if (params.caller == 'freebayes') {

        //VARIANT CALLING CON FREEBAYES

        //vcf  = Channel.empty()

        //FREEBAYES(ch_bam_mapped,fasta,fasta_fai)

        //vcf   = vcf.mix(FREEBAYES.out.vcf)

        //ch_versions   = ch_versions.mix(FREEBAYES.out.versions)


    //} else {

        //VARIANT CALLING CON STRELKA_GERMLINE

        //STRELKA_GERMLINE(ch_bam_mapped,ch_bam_bai,[],[], fasta, fasta_fai)

        //vcf = vcf.mix(STRELKA_GERMLINE.out.vcf)

        //ch_versions = ch_versions.mix(STRELKA_GERMLINE.out.versions)  

    //}

        //VARIANT CALLING CON MUTEC2
        //GATK4_MUTECT2(ch_bam_mapped,ch_bam_bai,intervals,fasta,fasta_fai,dict,germline_resource,germline_resource_tbi,[],[])


//*****************************************************************************************************************

    //BCFTOOLS

        //BCFTOOLS_STATS(vcf.map{meta, vcf -> [meta, vcf, []]}, [], [], [])

        //ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)


//*****************************************************************************************************************

    //VCFTOOLS ARREGLAR!!
        //No aparece ningun archivo en results
        //El bed puede no ser adecuado

        //VCFTOOLS(vcf,ch_bed_solo, [])

        //ch_versions = ch_versions.mix(VCFTOOLS.out.versions)

//*****************************************************************************************************************


       //RTGTOOLS_VCFEVAL()


//*****************************************************************************************************************

    //if (params.anotador == 'freebayes') {

        //ANOTACION CON SNPEFF

        //SNPEFF(vcf,snpeff_db,snpeff_cache)

        //ch_versions = ch_versions.mix(SNPEFF.out.versions)


    //} else {
    
        //ANOTACION CON VCF2MAF
        //Descomprimir vcf

        //UNCOMPRESS(vcf)

        //VCF2MAF(vcf,fasta,vep_cache)

    //}
        

//*****************************************************************************************************************


    //QC Opencga

        //OPENCGA_QC(ch_bam_mapped_bowtie2,ch_bam_bai,ch_bed,dict)


//*****************************************************************************************************************

    //ANOTACION CELLBASE

        //OPENCGA_ANNOTATE(freebayes_vcf)


//*****************************************************************************************************************


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml')
    )

//*****************************************************************************************************************

    // MULTIQC
    
    workflow_summary    = WorkflowAlberto.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowAlberto.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
