include wdl/.env
export

run_on_cori_jaws:
	jaws submit wdl/test.wdl storage/input_local.json TAHOMA

run_on_cori_local:
	java -Dconfig.file=cromwell_shifter.conf \
     -Dbackend.providers.Local.config.dockerRoot=$(pwd)/cromwell-executions \
     -jar /global/cfs/projectdirs/jaws/cromwell/cromwell.jar run wdl/test.wdl -i storage/input_local.json

run_test:
	java -jar wdl/cromwell/cromwell-66.jar run wdl/workflow/docker/main.wdl -i storage/input_local.json

validate_wdl:
	java -jar wdl/cromwell/womtool-66.jar validate wdl/workflow/docker/main.wdl -i storage/input_local.json

spit_dag:
	java -jar wdl/cromwell/womtool-66.jar womgraph wdl/workflow/docker/main.wdl

run_wdl:
	java -jar wdl/cromwell/cromwell-66.jar run wdl/workflow/docker/main.wdl -i storage/input_local.json

run_on_prismweb:
	~/cromwell/jdk-12.0.1/bin/java -jar ~/cromwell/cromwell-66.jar run wdl/workflow/docker/main.wdl -i storage/input.json

run_wdl_on_cori:
	java -jar /global/cfs/projectdirs/jaws/cromwell/cromwell-66.jar run wdl/workflow/shifter_at_cori/run_job_analysis.wdl -i /global/homes/a/anubhav/storage/input.json

prepare-your-input:
	wdl/scripts/prepare_input.sh

all-build:
	@docker-compose -f wdl/docker-compose.yml build

all-build-up:
	@docker-compose -f wdl/docker-compose.yml up --detach --build

just-pull-all:
	wdl/scripts/pull_spin_containers.sh

push-all-images:
	#@docker image tag wdl_msconvert:latest microbiomedata/metapro-msconvert:latest
	@docker push microbiomedata/metapro-msconvert:latest

	#@docker image tag wdl_masic_console:latest microbiomedata/metapro-masic:$(MASIC_VERSION)
	@docker push microbiomedata/metapro-masic:$(MASIC_VERSION)

	#@docker image tag wdl_mzidtotsvconverter:latest microbiomedata/metapro-mzidtotsvconverter:$(MZID2TSV_VERSION)
	@docker push microbiomedata/metapro-mzidtotsvconverter:$(MZID2TSV_VERSION)

	#@docker image tag wdl_msgfplus:latest microbiomedata/metapro-msgfplus:$(MSGFPLUS_VERSION)
	@docker push microbiomedata/metapro-msgfplus:$(MSGFPLUS_VERSION)

	#@docker image tag wdl_fastafilesplitter:latest microbiomedata/metapro-fastafilesplitter:$(Fasta_File_Splitter_VERSION)
	@docker push microbiomedata/metapro-fastafilesplitter:$(Fasta_File_Splitter_VERSION)

	#@docker image tag wdl_validatefastafile:latest microbiomedata/metapro-validatefastafile:$(VALIDATE_FASTA_FILE_VERSION)
	@docker push microbiomedata/metapro-validatefastafile:$(VALIDATE_FASTA_FILE_VERSION)

	#@docker image tag wdl_masicresultsmerge:latest microbiomedata/metapro-masicresultsmerge:$(MASICResultsMerge_VERSION)
	@docker push microbiomedata/metapro-masicresultsmerge:$(MASICResultsMerge_VERSION) 

	#@docker image tag wdl_proteindigestionsimulator:latest microbiomedata/metapro-proteindigestionsimulator:$(PDS_VERSION)
	@docker push microbiomedata/metapro-proteindigestionsimulator:$(PDS_VERSION)

	#@docker image tag wdl_mzidmerger:latest microbiomedata/metapro-mzidmerger:$(MzidMerger_VERSION)
	@docker push microbiomedata/metapro-mzidmerger:$(MzidMerger_VERSION)

	#@docker image tag wdl_peptidehitresultsprocrunner:latest microbiomedata/metapro-peptidehitresultsprocrunner:$(PHRP_VERSION)
	@docker push microbiomedata/metapro-peptidehitresultsprocrunner:$(PHRP_VERSION)

	#@docker image tag wdl_metadatacollection:latest microbiomedata/metapro-metadatacollection:2.0.0
	@docker push  microbiomedata/metapro-metadatacollection:2.0.0

	#@docker image tag wdl_post-processing:latest microbiomedata/metapro-post-processing:2.0.0
	@docker push microbiomedata/metapro-post-processing:2.0.0

#-------------
build_unified:
	docker-compose up --build
run_workflow:
	bash run_task.sh
#-------------