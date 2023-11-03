version 1.0

task convertToMgf {
    input {
        File   raw_file
        String mgf_file_name = basename(raw_file, ".raw")
    }
    command {
        wine msconvert \
            ${raw_file} \
            --zlib \
            --filter "peakPicking vendor msLevel=1-" \
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
        File kaiko_config
        String kaiko_volume_dir
    }
    command <<<
        execution_dir=$(pwd)
        cd /Kaiko_pipeline
        mgf=`basename ~{mgf_file}`
        mkdir input
        # ln -s ~{mgf_file} input/$mgf
        mv ~{mgf_file} input/$mgf
        ln -s ~{kaiko_volume_dir} Kaiko_volume
        python Kaiko_pipeline_main.py \
            --config=~{kaiko_config}
        find Kaiko_output -name '*.fasta' -exec mv -t $execution_dir {} +
        cd execution_dir
    >>>
    output {
        File outfile = glob("*.fasta")[0]
    }
    runtime {
        docker: 'kaiko-py310:latest'
    }
}

workflow run {
    String  kaiko_data_location="/refdata/Kaiko_volume/"
    input {
        File   raw_file
        File   kaiko_config
    }

    call convertToMgf {
        input:
            raw_file = raw_file,
    }
    call kaiko {
        input:
            mgf_file = convertToMgf.outfile,
            kaiko_config = kaiko_config,
            kaiko_volume_dir = kaiko_data_location
    }

    output {
        File   faa_file = kaiko.outfile
     }
}