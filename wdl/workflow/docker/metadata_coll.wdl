task collect{
    File resultant_file
    File peptide_file
    File protein_file
    File qc_metric_file
    String start_time
    String end_time
    String study

    command {
           python /app/metadata_collection/gen_metadata.py \
            ${resultant_file} \
            ${peptide_file} \
            ${protein_file} \
            ${qc_metric_file} \
            ${start_time} \
            ${end_time} \
            ${study}
    }
    output {
        File activity = "${study}_MetaProteomicAnalysis_activity.json"
        File data_object="${study}_emsl_analysis_data_objects.json"

    }
    runtime {
        docker: 'microbiomedata/metapro-metadatacollection:2.0.0'
    }
}
workflow gen_metadata{
    File resultant_file
    File peptide_file
    File protein_file
    File qc_metric_file
    String start_time
    String end_time
    String study

    call collect {
        input:
            resultant_file=resultant_file,
            peptide_file=peptide_file,
            protein_file=protein_file,
            qc_metric_file=qc_metric_file,
            start_time=start_time,
            end_time=end_time,
            study=study
    }
    output {
        File activity=collect.activity
        File data_object=collect.data_object
     }

}