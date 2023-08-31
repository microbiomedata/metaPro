version 1.0

task convertToMgf {
    input {
        File   raw_file
        String mgf_file_name = basename(raw_file, ".raw")
    }
    command {
        wine msconvert \
        ~{raw_file} \
        --mgf
    }
    output {
        File   outfile = "${mgf_file_name}.mgf"
    }
    runtime {
        docker: 'microbiomedata/metapro-msconvert:v3.0.21258'
    }
}

task kaiko {
    input {
        File mgf_file
    }
    command {
        mgf=`basename "${mgf_file}"`
        ln -s "${mgf_file}" "$mgf"
        echo '{"denovo": {"mgf_dir": "."}}' > overrides.json
        python run_kaiko_denovo.py \
            --config=overrides.json
    }
    output {
        File outfile = "${faa_file}.faa"
    }
    runtime {
        docker: 'kaiko-py310:latest'
    }
}

# task generateFeaturesAnnotations {
#     input {
#         File   faa_file
#     }
#     command {

#     }
#     output {
#         File   outfile = "${faa_file}.gff"
#     }
#     runtime {
#         docker: ''
#     }
# }

workflow report_gen{
    input{
        File   raw_file
        # String dataset_name
        File   kaiko_config
    }

    call convertToMgf {
        input:
            raw_file = raw_file,
            # dataset_name = dataset_name
    }
    call kaiko {
        input:
            raw_file = raw_file
            gff_file = gff_file,
    }
    # call generateFeaturesAnnotations {
    #     input:
    #         faa_file = kaiko.outfile
    # }

    output {
        File   faa_file = kaiko.outfile
        # File   gff_file = ficus_analysis.protein_file
     }

}