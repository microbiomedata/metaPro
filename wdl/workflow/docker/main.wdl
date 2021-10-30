import "run_job_analysis.wdl" as run_analysis
import "report_gen.wdl" as generate_reports
import "metadata_coll.wdl" as collect_metadata

workflow metapro {
    String study
    String dataset_id
    String genome_dir
    String dataset_name
    String annotation_name
    File raw_file_loc
    File faa_file_loc
    File gff_file_loc
    String QVALUE_THRESHOLD
    File MASIC_PARAM_FILENAME
    File MSGFPLUS_PARAM_FILENAME
    File CONTAMINANT_FILENAME

    call run_analysis.job_analysis {
        input:
            dataset_name=dataset_name,
            annotation_name=annotation_name,
            raw_file_loc=raw_file_loc,
            faa_file_loc=faa_file_loc,
            QVALUE_THRESHOLD=QVALUE_THRESHOLD,
            MASIC_PARAM_FILENAME=MASIC_PARAM_FILENAME,
            MSGFPLUS_PARAM_FILENAME=MSGFPLUS_PARAM_FILENAME,
            CONTAMINANT_FILENAME=CONTAMINANT_FILENAME
    }
    call generate_reports.report_gen {
        input:
            faa_txt_file=faa_file_loc,
            gff_file=gff_file_loc,
            resultant_file=job_analysis.resultant_file,
            Dataset_id=dataset_id,
            genome_directory=genome_dir,
            q_value_threshold=QVALUE_THRESHOLD,
            annotation_name=annotation_name

    }

    call collect_metadata.gen_metadata {
        input:
            study=study,
            resultant_file=job_analysis.resultant_file,
            peptide_file=report_gen.peptide_file,
            protein_file=report_gen.protein_file,
            qc_metric_file=report_gen.qc_metric_file,
            start_time=job_analysis.start_time,
            end_time=job_analysis.end_time
    }
}
