version 1.0

task masic {
    input{
        File   raw_file
        File   masic_param
        String dataset_name
    }
    command {
        mono /app/masic/MASIC_Console.exe \
        /I:~{raw_file} \
        /O:~{'.'} \
        /P:~{masic_param}
    }
    output {
        File   outfile = "${dataset_name}_SICstats.txt"
    }
    runtime {
        docker: 'microbiomedata/metapro-masic:v3.2.7901'
    }
}
task msconvert {
    input{
        File   raw_file
        String dataset_name
    }
    command {
        wine msconvert \
        ~{raw_file} \
        --zlib \
        --filter 'peakPicking true 1-'
    }
    output {
        File   outfile = "${dataset_name}.mzML"
    }
    runtime {
        docker: 'microbiomedata/metapro-msconvert:v3.0.21258'
    }
}
task msgfplus {
    input{
        File   mzml_file
        File   contaminated_fasta_file
        File   msgfplus_params
        String dataset_name
        String revcat_name = basename(sub(contaminated_fasta_file, "\\.fasta$", ".revCat.fasta"))
    }
    command {
        java -Xmx32G -jar /app/msgf/MSGFPlus.jar \
        -s ~{mzml_file} \
        -d ~{contaminated_fasta_file} \
        -o "~{dataset_name}.mzid" \
        -conf ~{msgfplus_params} \
        -thread 16 \
        -verbose 1
        echo '>>> moving revCat.fasta file in execution folder.'
        rev_cat_fasta_loc=$(find .. -type f -regex ".*~{revcat_name}")
        cp $rev_cat_fasta_loc ../execution/
    }
    output {
        File   outfile       = "${dataset_name}.mzid"
        File   rev_cat_fasta = revcat_name
    }
    runtime {
        docker: 'microbiomedata/metapro-msgfplus:v2022.04.18'
    }
}
task mzidtotsvconverter{
    input{
        File   mzid_file
        String dataset_name
    }
    command {
        mono /app/mzid2tsv/net462/MzidToTsvConverter.exe \
        -mzid:~{mzid_file} \
        -tsv:"~{dataset_name}.tsv" \
        -unroll \
        -showDecoy
        echo '>>> moving tsv file in execution folder.'
        tsv_file_loc=$(find .. -type f -regex ".*~{dataset_name}.tsv")
        cp $tsv_file_loc ../execution/
    }
    output {
        File   outfile = "${dataset_name}.tsv"
    }
    runtime {
        docker: 'microbiomedata/metapro-mzidtotsvconverter:v1.4.6'
    }
}
task peptidehitresultsprocrunner {
    input{
        File   tsv_file
#        File   msgfplus_modef_params
#        File   mass_correction_params
        File   msgfplus_params
        File   revcatfasta_file
        String dataset_name
    }
#        -M:~{msgfplus_modef_params} \
#        -T:~{mass_correction_params} \
    command {
        mono /app/phrp/PeptideHitResultsProcRunner.exe \
        -I:~{tsv_file} \
        -N:~{msgfplus_params} \
        -SynPvalue:0.2 \
        -SynProb:0.05 \
        -ProteinMods \
        -F:~{revcatfasta_file} \
        -O:~{'.'}
    }
    output {
        File   outfile = "${dataset_name}_syn.txt"
        File   first_hits_file = "${dataset_name}_fht.txt"
    }
    runtime {
        docker: 'microbiomedata/metapro-peptidehitresultsprocrunner:v3.0.7842'
    }
}
task masicresultmerge {
    input{
        File   sic_stats_file
        File   synopsis_file
        String dataset_name
        String dataset_id
        String faa_file_id
    }
    command {
        synopsis_file_loc=$(find .. -type f -regex ".*~{dataset_name}_syn.txt")
        cp $synopsis_file_loc ../execution/
        sic_stats_file_loc=$(find .. -type f -regex ".*~{dataset_name}_SICstats.txt")
        cp $sic_stats_file_loc ../execution/
        mv ~{dataset_name}_syn.txt ~{dataset_id}_~{faa_file_id}_msgfplus_syn.txt
        mv ~{dataset_name}_SICstats.txt ~{dataset_id}_~{faa_file_id}_SICStats.txt
        mono /app/MASICResultsMerge/MASICResultsMerger.exe \
        ~{dataset_id}_~{faa_file_id}_msgfplus_syn.txt
        date --iso-8601=seconds > stop.txt
    }
    output {
        File   outfile = "${dataset_id}_${faa_file_id}_msgfplus_syn_PlusSICStats.txt"
        String stop = read_string("stop.txt")
    }
    runtime {
        docker: 'microbiomedata/metapro-masicresultsmerge:v2.0.7983'
    }
}
task fastaFileSplitter {
    input{
        File    fasta_file_loc
        Int     split_size
    }
    command {
        mono /app/FastaFileSplitter/FastaFileSplitter.exe \
            /I:~{fasta_file_loc} \
            /MB:~{split_size} \
            /O:~{'.'}
    }
    output {
        Array[File]   outfiles = glob("*.fasta")
    }
    runtime {
        docker: 'microbiomedata/metapro-fastafilesplitter:v1.1.7887'
    }
}
task msgfplusresultsmerge {
    input{
        Array[File] mzid_files
        Array[File] fasta_files
        String output_mzid_file_name
        String output_fasta_file_name
    }
    command<<<
        fps=( ~{sep=' ' mzid_files } )
        for i in "${!fps[@]}"
        do
            ln ${fps[$i]} `basename ${fps[$i]}`
        done
        mono /app/MzidMerger/net472/MzidMerger.exe \
            -inDir:~{'.'} \
            -filter:"*.mzid" \
            -out:~{output_mzid_file_name} \
            -keepOnlyBestResults
        cat ~{sep=' ' fasta_files } >> ~{output_fasta_file_name}
    >>>
    output {
        File    outfile_mzid = output_mzid_file_name
        File    outfile_fasta = output_fasta_file_name 
    }
    runtime {
        docker: 'microbiomedata/metapro-mzidmerger:v1.4.26'
    }
}
task concatcontaminate {
    input{
        File faa_file
        File contaminate_file
        String output_filename = basename(faa_file)
    }
    command<<<
        date -u --iso-8601=seconds > start.txt
        cat ~{faa_file} ~{contaminate_file} > ~{output_filename}
    >>>
    output {
        File    outfile = output_filename
        String  start = read_string("start.txt")
    }
    runtime {
        docker: 'microbiomedata/metapro-masic:v3.2.7901'
    }
}

workflow job_analysis{
    input{
        String dataset_name
        File   raw_file_loc
        File   faa_file_loc
        String QVALUE_THRESHOLD
        File   MASIC_PARAM_FILENAME
        File   MSGFPLUS_PARAM_FILENAME
        File   CONTAMINANT_FILENAME
        Int    FASTA_SPLIT_ON_SIZE_MB
        String dataset_id
        String faa_file_id
    }
    call concatcontaminate {
        input:
            faa_file        = faa_file_loc,
            contaminate_file= CONTAMINANT_FILENAME
    }
    call masic {
        input:
            raw_file    = raw_file_loc,
            masic_param = MASIC_PARAM_FILENAME,
            dataset_name= dataset_name
    }
    call msconvert {
        input:
            raw_file     = raw_file_loc,
            dataset_name = dataset_name
    }
    if(size(faa_file_loc, 'MB') > FASTA_SPLIT_ON_SIZE_MB)
    {
        call fastaFileSplitter {
            input:
                fasta_file_loc  = faa_file_loc,
                split_size  = FASTA_SPLIT_ON_SIZE_MB
        }
        scatter(split_fasta_file in fastaFileSplitter.outfiles)
        {
            call concatcontaminate as concatcontaminatesplit {
                input:
                    faa_file        = split_fasta_file,
                    contaminate_file= CONTAMINANT_FILENAME
            }
            # Number output MS-GF+ .mzid files similarly to input FASTAs for recordkeeping
            String faa_file_basename = basename(faa_file_loc, ".faa")
            String numbered_dataset_name = sub(split_fasta_file, faa_file_basename, dataset_name)
            String dataset_basename = basename(numbered_dataset_name, ".fasta")

            call msgfplus as msgfplussplit{
                input:
                    mzml_file               = msconvert.outfile,
                    contaminated_fasta_file = concatcontaminatesplit.outfile,
                    msgfplus_params         = MSGFPLUS_PARAM_FILENAME,
                    dataset_name            = dataset_basename,
            }
        }
        call msgfplusresultsmerge {
            input:
                mzid_files = msgfplussplit.outfile,
                fasta_files = msgfplussplit.rev_cat_fasta,
                output_mzid_file_name = "~{dataset_name}.mzid",
                output_fasta_file_name = basename(faa_file_loc, ".faa") + ".revCat.fasta"
        }

        File? msgfplus_split_and_merged = msgfplusresultsmerge.outfile_mzid
        File? rev_cat_fasta_split_and_merged = msgfplusresultsmerge.outfile_fasta
        Boolean? fasta_size_greater_than_split_size = true
    }
    if(size(faa_file_loc, 'MB') <= FASTA_SPLIT_ON_SIZE_MB)
    {
        call msgfplus{
            input:
                mzml_file               = msconvert.outfile,
                contaminated_fasta_file = concatcontaminate.outfile,
                msgfplus_params         = MSGFPLUS_PARAM_FILENAME,
                dataset_name            = dataset_name,
        }

        File? msgfplus_not_split = msgfplus.outfile
        File? rev_cat_fasta_not_split = msgfplus.rev_cat_fasta
        Boolean? fasta_size_less_than_or_equal_split_size = false
    }

    # These arrays should have only a single valid value
    Array[File?] msgfplus_mzids          = [msgfplus_not_split, msgfplus_split_and_merged]
    Array[File?] msgfplus_revCata_fastas = [rev_cat_fasta_not_split, rev_cat_fasta_split_and_merged]

    # The value of the single valid boolean indicates if we've taken the split FASTA processing route
    Array[Boolean?] fasta_split_state    = [fasta_size_greater_than_split_size, fasta_size_less_than_or_equal_split_size]

    call mzidtotsvconverter {
        input:
            mzid_file    = select_first(msgfplus_mzids),
            dataset_name = dataset_name
    }
    call peptidehitresultsprocrunner {
        input:
            tsv_file               = mzidtotsvconverter.outfile,
#            msgfplus_modef_params  = "",
#            mass_correction_params = "",
            msgfplus_params        = MSGFPLUS_PARAM_FILENAME,
            revcatfasta_file       = select_first(msgfplus_revCata_fastas),
            dataset_name           = dataset_name
    }
    call masicresultmerge {
        input:
            sic_stats_file = masic.outfile,
            synopsis_file  = peptidehitresultsprocrunner.outfile,
            dataset_name   = dataset_name,
            dataset_id     = dataset_id,
            faa_file_id    = faa_file_id
    }

    output {
        File   resultant_file = masicresultmerge.outfile
        File   faa_with_contaminates = concatcontaminate.outfile
        File   first_hits_file = peptidehitresultsprocrunner.first_hits_file
        String start_time= concatcontaminate.start
        String end_time=masicresultmerge.stop
        Boolean did_split = select_first(fasta_split_state)
     }
}
