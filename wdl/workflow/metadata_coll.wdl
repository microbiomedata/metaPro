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
        String pipeline_type
        String execution_resource
        String git_url
    }
    command {
        python /app/metadata_collection/code/gen_metadata.py \
            ~{write_json(results)} \
            "~{study}" \
            ~{pipeline_type} \
            "~{execution_resource}" \
            ~{git_url} \
            "nmdc:DataObject"
    }
    output {
        File   activity    = "${study}_MetaProteomicAnalysis_activity.json"
        File   data_object = "${study}_emsl_analysis_data_objects.json"
    }
    runtime {
        docker: 'microbiomedata/metapro-metadatacollection:2.1.0'
    }
}
workflow gen_metadata{
    input{
        String study
        Array[Result] results
        String pipeline_type
        String execution_resource
        String git_url
    }
    call collect {
        input:
            study          = study,
            results        = results,
            pipeline_type  = pipeline_type, 
            execution_resource = execution_resource,
            git_url = git_url
    }
    output {
        File   activity    = collect.activity
        File   data_object = collect.data_object
     }
}