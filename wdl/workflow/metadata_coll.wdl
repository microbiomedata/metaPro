version 1.0

struct Result {
    File resultant_file
    File peptide_report_file
    File protein_report_file
    File qc_metric_report_file
    File faa_file
    File txt_faa_file
    File gff_file
    String genome_directory
    String dataset_id
    String start_time
    String end_time
    String fasta_id
    String gff_id
}

task collect{
    input{
        String study
        Array[Result] results
        String execution_resource
        String git_url
        String data_url
        File contaminate_file
        File masic_param_file
        File msgfplus_param_file
        String contaminant_file_id
        String masic_param_id
        String msgfplus_param_id
        String version
        Boolean metagenome_free
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
            ~{msgfplus_param_file} \
            ~{masic_param_id} \
            ~{msgfplus_param_id} \
            ~{contaminant_file_id} \
            ~{version} \
            ~{metagenome_free}
    }
    output {
        File   activity    = "${study}_MetaProteomicAnalysis_activity.json"
        File   data_object = "${study}_analysis_data_objects.json"
    }
    runtime {
        docker: 'microbiomedata/metapro-metadatacollection:1.2.1'
    }
}
workflow gen_metadata{
    input{
        String study
        Array[Result] results
        String execution_resource
        String git_url
        String data_url
        File contaminate_file
        File masic_param_file
        File msgfplus_param_file
        String contaminant_file_id
        String masic_param_id
        String msgfplus_param_id
        String version
        Boolean metagenome_free
    }
    call collect {
        input:
            study = study,
            results = results,
            execution_resource = execution_resource,
            git_url = git_url,
            data_url = data_url,
            contaminate_file = contaminate_file,
            masic_param_file = masic_param_file,
            msgfplus_param_file = msgfplus_param_file,
            contaminant_file_id = contaminant_file_id,
            masic_param_id = masic_param_id,
            msgfplus_param_id = msgfplus_param_id,
            version = version,
            metagenome_free = metagenome_free
    }
    output {
        File   activity    = collect.activity
        File   data_object = collect.data_object
     }
}