version 1.0

struct Result {
    File resultant_file
    File peptide_report_file
    File protein_report_file
    File qc_metric_report_file
    File faa_file
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
        String data_url
        File contaminate_file
        File masic_param_file
        File msgfplus_param_file
    }
    command {
        python /app/metadata_collection/code/gen_metadata.py \
            ~{write_json(results)} \
            "~{study}" \
            "~{execution_resource}" \
            ~{git_url} \
            ~{data_url} \
            ~{contaminate_file} \
            ~{masic_param_file} \
            ~{msgfplus_param_file}
    }
    output {
        File   activity    = "${study}_MetaProteomicAnalysis_activity.json"
        File   data_object = "${study}_analysis_data_objects.json"
    }
    runtime {
        docker: 'microbiomedata/metapro-metadatacollection:1.1.0'
    }
}
workflow gen_metadata{
    input{
        String study
        Array[Result] results
        String pipeline_type
        String execution_resource
        String git_url
        String data_url
        File contaminate_file
        File masic_param_file
        File msgfplus_param_file
    }
    call collect {
        input:
            study = study,
            results = results,
            pipeline_type = pipeline_type, 
            execution_resource = execution_resource,
            git_url = git_url,
            data_url = data_url,
            contaminate_file = contaminate_file,
            masic_param_file = masic_param_file,
            msgfplus_param_file = msgfplus_param_file
    }
    output {
        File   activity    = collect.activity
        File   data_object = collect.data_object
     }
}