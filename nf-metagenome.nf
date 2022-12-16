#!/usr/bin/env nextflow

/* enable DSL2*/
nextflow.enable.dsl=2

/*

/*
================================================================
The AGUILAR LAB presents...

  The Metagenome Assembly Pipeline


==================================================================
Version: 0.1
Project repository:
==================================================================
Authors:

- Bioinformatics Design


- Bioinformatics Development


- Nextflow Port
 Israel Aguilar-Ordonez (iaguilaror@gmail.com)

=============================
Pipeline Processes In Brief:
=============================

Pre-processing:


Core-processing:


Pos-processing


Anlysis


ENDING
 _register_configs

================================================================*/

/* Define the help message as a function to call when needed *//////////////////////////////
def helpMessage() {
	log.info"""
  ==========================================
							${pipeline_name}

	The global aligment pipeline
  v${version}
  ==========================================

	Usage:

	nextflow run ${pipeline_name}.nf --fastq_dir <path to input 1>  [--output_dir path to results ]

		--ref_fa	<-  Directory whith compressed FASTQ files;

	  --output_dir     <- directory where results, intermediate and log files will be stored;
	      default: same dir where --query_fasta resides

	  -resume	   <- Use cached results if the executed project has been run before;
	      default: not activated
	      This native NF option checks if anything has changed from a previous pipeline execution.
	      Then, it resumes the run from the last successful stage.
	      i.e. If for some reason your previous run got interrupted,
	      running the -resume option will take it from the last successful pipeline stage
	      instead of starting over
	      Read more here: https://www.nextflow.io/docs/latest/getstarted.html#getstart-resume
	  --help           <- Shows Pipeline Information
	  --version        <- Show version
	""".stripIndent()
}

/*//////////////////////////////
  Define pipeline version
  If you bump the number, remember to bump it in the header description at the begining of this script too
*/
version = "0.1"

/*//////////////////////////////
  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
pipeline_name = "nf-metagenome"

/*
  Initiate default values for parameters
  to avoid "WARN: Access to undefined parameter" messages
*/
params.fastq_dir = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.help = false //default is false to not trigger help message automatically at every run
params.version = false //default is false to not trigger version message automatically at every run

/*//////////////////////////////
  If the user inputs the --help flag
  print the help message and exit pipeline
*/
if (params.help){
	helpMessage()
	exit 0
}

/*//////////////////////////////
  If the user inputs the --version flag
  print the pipeline version
*/
if (params.version){
	println "${pipeline_name} v${version}"
	exit 0
}

/*//////////////////////////////
  Define the Nextflow version under which this pipeline was developed or successfuly tested
  Updated by iaguilar at MAY 2021
*/
nextflow_required_version = '20.01.0'
/*
  Try Catch to verify compatible Nextflow version
  If user Nextflow version is lower than the required version pipeline will continue
  but a message is printed to tell the user maybe it's a good idea to update her/his Nextflow
*/
try {
	if( ! nextflow.version.matches(">= $nextflow_required_version") ){
		throw GroovyException('Your Nextflow version is older than Pipeline required version')
	}
} catch (all) {
	log.error "-----\n" +
			"  This pipeline requires Nextflow version: $nextflow_required_version \n" +
      "  But you are running version: $workflow.nextflow.version \n" +
			"  The pipeline will continue but some things may not work as intended\n" +
			"  You may want to run `nextflow self-update` to update Nextflow\n" +
			"============================================================"
}

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

/* Check if the input directory is provided
    if it was not provided, it keeps the 'false' value assigned in the parameter initiation block above
    and this test fails
*/
if ( !params.fastq_dir) {
  log.error " Please provide the --fastq_dir \n\n" +
  " For more information, execute: nextflow run ${pipeline_name}.nf --help"
  exit 1
}

/*
Output directory definition
Default value to create directory is the parent dir of --input_dir
*/
params.output_dir = file(params.fastq_dir).getParent() //!! maybe creates bug, should check

/*
  Results and Intermediate directory definition
  They are always relative to the base Output Directory
  and they always include the pipeline name in the variable pipeline_name defined by this Script

  This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*
========================================================================================
    LOG RUN INFORMATION
========================================================================================
*/
log.info"""
==========================================
The ${pipeline_name} pipeline
v${version}
==========================================
"""
log.info "--Nextflow metadata--"
/* define function to store nextflow metadata summary info */
def nfsummary = [:]
/* log parameter values beign used into summary */
/* For the following runtime metadata origins, see https://www.nextflow.io/docs/latest/metadata.html */
nfsummary['Resumed run?'] = workflow.resume
nfsummary['Run Name']			= workflow.runName
nfsummary['Current user']		= workflow.userName
/* string transform the time and date of run start; remove : chars and replace spaces by underscores */
nfsummary['Start time']			= workflow.start.toString().replace(":", "").replace(" ", "_")
nfsummary['Script dir']		 = workflow.projectDir
nfsummary['Working dir']		 = workflow.workDir
nfsummary['Current dir']		= workflow.launchDir
nfsummary['Launch command'] = workflow.commandLine
log.info nfsummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "\n\n--Pipeline Parameters--"
/* define function to store nextflow metadata summary info */
def pipelinesummary = [:]
/* log parameter values beign used into summary */
pipelinesummary['input directory']						= params.fastq_dir
pipelinesummary['Results Dir']								= results_dir
pipelinesummary['Intermediate Dir']						= intermediates_dir
/* print stored summary info */
log.info pipelinesummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================\nPipeline Start"

////////////////////////////////////////////////////
/* --  READ INPUTS  -- */
////////////////////////////////////////////////////

/* Define function for finding files that share sample name */
/* in this case, the file name comes from the 1st ".",
since tokenize array starts at 0, array index shoould be 0 */
def get_sample_prefix = { file -> file.name}

/* Load fq files into channel as pairs */
Channel
  .fromPath( "${params.fastq_dir}" )
  .set{fasta_inputs}

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
									/* PRE-processing */

																		/* pos-processing */
												/*======== GTDBTK  ========*/
include {GTDBTK as GTDBTK_METASPADES } from './modules/local/pos/gtdbtk/main.nf' addParams(tool: "metaspades")
include {GTDBTK as GTDBTK_MEGAHIT } from './modules/local/pos/gtdbtk/main.nf' addParams(tool: "megahit")


/*
========================================================================================
    											RUN MAIN WORKFLOW
========================================================================================
*/
workflow  {
					FASTA_BINS_MEGAHIT = fasta_inputs
				// GTDBTK
				GTDBTK_MEGAHIT(FASTA_BINS_MEGAHIT)
}
