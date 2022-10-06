version 1.0

struct Result {
    File resultant_file
    File peptide_report_file
    File protein_report_file
    File qc_metric_report_file
    File faa_file
    File contaminate_file
    File txt_faa_file
    String genome_directory
    String dataset_id
    String start_time
    String end_time
}

task collect{
    input{
        String study
        Array[Result] results
    }
    command {
        python /app/metadata_collection/code/gen_metadata.py \
            ${write_json(results)}
    }
    output {
        File   activity    = "${study}_MetaProteomicAnalysis_activity.json"
        File   data_object = "${study}_emsl_analysis_data_objects.json"
    }
    runtime {
        docker: 'microbiomedata/metapro-metadatacollection:2.0.1'
    }
}
workflow gen_metadata{
    input{
        String study
        Array[Result] results
    }
    call collect {
        input:
            study          = study,
            results        = results
    }
    output {
        File   activity    = collect.activity
        File   data_object = collect.data_object
     }
}