version 1.0

task ficus_analysis {
    input{
        File   faa_txt_file
        File   gff_file
        File   resultant_file
        File   first_hits_file
        String Dataset_id
        String faa_file_id
        String q_value_threshold
        String dataset_name
        Boolean did_split
        Boolean metagenome_free
    }
    command {
        if [ ${did_split} == true ]; then
            python /app/post-processing/compute_fdr.py \
                --file=${first_hits_file}   \
                --fdr=${q_value_threshold}  \
                --spece=1.0e-10 \
                --speceinc=1.0e-13 
            threshold=$(cat out.json | jq ".SpecEValue")
        else
            threshold=${q_value_threshold}
        fi
        python /app/post-processing/ficus_analysis.py \
            ${faa_txt_file} \
            ${gff_file} \
            ${resultant_file} \
            ${Dataset_id} \
            ${faa_file_id} \
            $threshold \
            ${dataset_name} \
            ${did_split} \
            ${metagenome_free}
    }
    output {
        File   peptide_file   = "${Dataset_id}_${faa_file_id}_Peptide_Report.tsv"
        File   protein_file   = "${Dataset_id}_${faa_file_id}_Protein_Report.tsv"
        File   qc_metric_file = "${Dataset_id}_${faa_file_id}_QC_metrics.tsv"
    }
    runtime {
        docker: 'microbiomedata/metapro-post-processing:2.0.0'
    }
}
task proteinDigestionSimulator {
    input{
        File faa_file
        String output_filename = sub(basename(faa_file), "\\.faa$", ".txt")
    }
    command {
        mono /app/pds/ProteinDigestionSimulator.exe \
        -I:~{faa_file} \
        -F \
        -O:~{'.'}
    }
    output {
        File outfile = output_filename
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
        File   first_hits_file
        String Dataset_id
        String faa_file_id
        String q_value_threshold
        String dataset_name
        Boolean did_split
        Boolean metagenome_free
    }

    call proteinDigestionSimulator {
        input:
            faa_file        = faa_txt_file,
    }
    
    call ficus_analysis {
        input:
            faa_txt_file      = proteinDigestionSimulator.outfile,
            gff_file          = gff_file,
            resultant_file    = resultant_file,
            first_hits_file   = first_hits_file,
            Dataset_id        = Dataset_id,
            faa_file_id       = faa_file_id,
            q_value_threshold = q_value_threshold,
            dataset_name      = dataset_name,
            did_split         = did_split,
            metagenome_free   = metagenome_free
    }
    output {
        File   peptide_file   = ficus_analysis.peptide_file
        File   protein_file   = ficus_analysis.protein_file
        File   qc_metric_file = ficus_analysis.qc_metric_file
        File   txt_faa_file   = proteinDigestionSimulator.outfile
     }

}
