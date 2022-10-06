version 1.0

import "run_job_analysis.wdl" as run_analysis
import "report_gen.wdl" as generate_reports
import "metadata_coll.wdl" as collect_metadata

workflow metapro {

    input{
        Array[Object] mapper_list
        String QVALUE_THRESHOLD
        File MASIC_PARAM_FILE_LOC
        File MSGFPLUS_PARAM_FILE_LOC
        File CONTAMINANT_FILE_LOC
        String STUDY
        Int FASTA_SPLIT_ON_SIZE_MB
        Int FASTA_SPLIT_COUNT
    }

    scatter (myobj in mapper_list) {
        call run_analysis.job_analysis {
            input:
                dataset_name            = myobj['dataset_name'],
                annotation_name         = myobj['annotation_name'],
                raw_file_loc            = myobj['raw_file_loc'],
                faa_file_loc            = myobj['faa_file_loc'],
                QVALUE_THRESHOLD        = QVALUE_THRESHOLD,
                MASIC_PARAM_FILENAME    = MASIC_PARAM_FILE_LOC,
                MSGFPLUS_PARAM_FILENAME = MSGFPLUS_PARAM_FILE_LOC,
                CONTAMINANT_FILENAME    = CONTAMINANT_FILE_LOC,
                FASTA_SPLIT_ON_SIZE_MB  = FASTA_SPLIT_ON_SIZE_MB,
                FASTA_SPLIT_COUNT       = FASTA_SPLIT_COUNT
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
                dataset_name      = myobj['dataset_name']
        }

        Result result = {
            "resultant_file": job_analysis.resultant_file,
            "peptide_report_file": report_gen.peptide_file,
            "protein_report_file": report_gen.protein_file,
            "qc_metric_report_file": report_gen.qc_metric_file,
            "faa_file": myobj['faa_file_loc'],
            "contaminate_file": CONTAMINANT_FILE_LOC,
            "txt_faa_file": report_gen.txt_faa_file,
            "genome_directory": myobj['genome_dir'],
            "dataset_id": myobj['dataset_id'],
            "start_time": job_analysis.start_time,
            "end_time": job_analysis.end_time
        }
    }

    Array[Result?] results_maybe = result
    Array[Result] results = select_all(results_maybe)

    call collect_metadata.gen_metadata {
        input:
            study       = STUDY,
            results     = results
    }
}
