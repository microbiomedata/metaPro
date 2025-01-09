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
        mkdir -p /data/input
        mkdir /data/output
        mgf=`basename ~{mgf_file}`
        ln -s ~{mgf_file} /data/input/$mgf
        cd /Kaiko_metaproteome
        python -m Kaiko_main.Kaiko_main \
            --config=~{kaiko_config}
        find /data/output/Kaiko_output -name '*.fasta' -exec mv -t $execution_dir {} +
        find /data/output/Kaiko_output -name '*.gff' -exec mv -t $execution_dir {} +
        cd $execution_dir
    >>>
    output {
        File outfile_fasta = glob("*.fasta")[0]
        File outfile_gff = glob("*.gff")[0]
    }
    runtime {
        docker: 'camiloposso15/kaiko_2.0-py3.10:latest'
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
        File   faa_file = kaiko.outfile_fasta
        File   gff_file = kaiko.outfile_gff
     }
}