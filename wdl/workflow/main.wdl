version 1.0

import "run_job_analysis.wdl" as run_analysis
import "report_gen.wdl" as generate_reports
import "metadata_coll.wdl" as collect_metadata

workflow metapro {
    Int fasta_split_on_size_mb = 1000
    Int fasta_split_count = 22
    String git_url = "https://github.com/microbiomedata/metaPro/releases/tag/2.0.0"

    input{
        Array[Object] mapper_list
        String QVALUE_THRESHOLD
        File MASIC_PARAM_FILE_LOC
        String MASIC_PARAM_FILE_ID
        File MSGFPLUS_PARAM_FILE_LOC
        String MSGFPLUS_PARAM_FILE_ID
        File CONTAMINANT_FILE_LOC
        String CONTAMINANT_FILE_ID
        String STUDY
        String EXECUTION_RESOURCE
        String DATA_URL
    }

    scatter (myobj in mapper_list) {
        call run_analysis.job_analysis {
            input:
                dataset_name            = myobj['dataset_name'],
                raw_file_loc            = myobj['raw_file_loc'],
                faa_file_loc            = myobj['faa_file_loc'],
                QVALUE_THRESHOLD        = QVALUE_THRESHOLD,
                MASIC_PARAM_FILENAME    = MASIC_PARAM_FILE_LOC,
                MSGFPLUS_PARAM_FILENAME = MSGFPLUS_PARAM_FILE_LOC,
                CONTAMINANT_FILENAME    = CONTAMINANT_FILE_LOC,
                FASTA_SPLIT_ON_SIZE_MB  = fasta_split_on_size_mb,
                FASTA_SPLIT_COUNT       = fasta_split_count
        }
        call generate_reports.report_gen {
            input:
                faa_txt_file      = job_analysis.faa_with_contaminates,
                gff_file          = myobj['gff_file_loc'],
                resultant_file    = job_analysis.resultant_file,
                Dataset_id        = myobj['dataset_id'],
                genome_directory  = myobj['genome_dir'],
                q_value_threshold = QVALUE_THRESHOLD,
                annotation_name   = myobj['annotation_name'],
                dataset_name      = myobj['dataset_name'],
                first_hits_file   = job_analysis.first_hits_file,
                did_split         = job_analysis.did_split
        }

        Result result = {
            "resultant_file": job_analysis.resultant_file,
            "peptide_report_file": report_gen.peptide_file,
            "protein_report_file": report_gen.protein_file,
            "qc_metric_report_file": report_gen.qc_metric_file,
            "faa_file": myobj['faa_file_loc'],
            "txt_faa_file": report_gen.txt_faa_file,
            "genome_directory": myobj['genome_dir'],
            "dataset_id": myobj['dataset_id'],
            "start_time": job_analysis.start_time,
            "end_time": job_analysis.end_time,
            "fasta_id": myobj['faa_file_id'],
            "gff_id": myobj['gff_file_id']
        }
    }

    Array[Result?] results_maybe = result
    Array[Result] results = select_all(results_maybe)

    call collect_metadata.gen_metadata {
        input:
            study       = STUDY,
            results     = results,
            execution_resource = EXECUTION_RESOURCE,
            git_url = git_url,
            data_url = DATA_URL,
            contaminate_file = CONTAMINANT_FILE_LOC,
            masic_param_file = MASIC_PARAM_FILE_LOC,
            msgfplus_param_file = MSGFPLUS_PARAM_FILE_LOC,
            contaminant_file_id = CONTAMINANT_FILE_ID,
            masic_param_id = MASIC_PARAM_FILE_ID,
            msgfplus_param_id = MSGFPLUS_PARAM_FILE_ID
    }
}
