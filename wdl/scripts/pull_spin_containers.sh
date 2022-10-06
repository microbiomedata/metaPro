#!/bin/bash
export $(grep -v '^#' wdl/.env | xargs)

images=(microbiomedata/metapro-msconvert:latest \
        microbiomedata/metapro-masic:$MASIC_VERSION \
        microbiomedata/metapro-mzidtotsvconverter:$MZID2TSV_VERSION \
        microbiomedata/metapro-msgfplus:$MSGFPLUS_VERSION \
        microbiomedata/metapro-peptidehitresultsprocrunner:$PHRP_VERSION \
        microbiomedata/metapro-mzidmerger:$MzidMerger_VERSION \
        microbiomedata/metapro-proteindigestionsimulator:$PDS_VERSION \
        microbiomedata/metapro-masicresultsmerge:$MASICResultsMerge_VERSION \
        microbiomedata/metapro-validatefastafile:$VALIDATE_FASTA_FILE_VERSION \
        microbiomedata/metapro-fastafilesplitter:$Fasta_File_Splitter_VERSION \
        microbiomedata/metapro-post-processing:2.0.0 \
        microbiomedata/metapro-metadatacollection:2.0.1 )

for image in ${images[*]} ; do
    docker pull $image;
    docker run -v /$(pwd)/storage:/app/storage:rw \
            -v /$(pwd)/logs:/app/logs:rw \
            $image;
  done

#  --detach -ti \
#            --name $(awk '{ sub(/.*\//, ""); sub(/:.*/, ""); print }' <<< "$image") \