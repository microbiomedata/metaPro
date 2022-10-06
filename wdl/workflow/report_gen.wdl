version 1.0

task ficus_analysis {
    input{
        File   faa_txt_file
        File   gff_file
        File   resultant_file
        String Dataset_id
        String genome_directory
        String q_value_threshold
        String dataset_name
    }
    command {
       python /app/post-processing/ficus_analysis.py \
        ~{faa_txt_file} \
        ~{gff_file} \
        ~{resultant_file} \
        ~{Dataset_id} \
        ~{genome_directory} \
        ~{q_value_threshold} \
        ~{dataset_name}
    }
    output {
        File   peptide_file   = "${Dataset_id}_${genome_directory}_Peptide_Report.tsv"
        File   protein_file   = "${Dataset_id}_${genome_directory}_Protein_Report.tsv"
        File   qc_metric_file = "${Dataset_id}_${genome_directory}_QC_metrics.tsv"
    }
    runtime {
        docker: 'microbiomedata/metapro-post-processing:2.0.0'
    }
}
task proteinDigestionSimulator {
    input{
        File faa_file
        String annotation_name
    }
    command {
        mono /app/pds/ProteinDigestionSimulator.exe \
        -I:~{faa_file} \
        -F \
        -O:~{'.'}
    }
    output {
        File outfile = "${annotation_name}_proteins.txt"
    }
    runtime {
        docker: "microbiomedata/metapro-proteindigestionsimulator:v2.3.7794"
    }
}

workflow report_gen{
    input{
        File   faa_txt_file
        File   gff_file
        File   resultant_file
        String Dataset_id
        String genome_directory
        String annotation_name
        String q_value_threshold
        String dataset_name
    }

    call proteinDigestionSimulator {
        input:
            faa_file        = faa_txt_file,
            annotation_name = annotation_name
    }
    call ficus_analysis {
        input:
            faa_txt_file      = proteinDigestionSimulator.outfile,
            gff_file          = gff_file,
            resultant_file    = resultant_file,
            Dataset_id        = Dataset_id,
            genome_directory  = genome_directory,
            q_value_threshold = q_value_threshold,
            dataset_name      = dataset_name
    }
    output {
        File   peptide_file   = ficus_analysis.peptide_file
        File   protein_file   = ficus_analysis.protein_file
        File   qc_metric_file = ficus_analysis.qc_metric_file
        File   txt_faa_file   = proteinDigestionSimulator.outfile
     }

}