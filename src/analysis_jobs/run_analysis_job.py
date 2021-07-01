import os
import fnmatch
import json
from utility.utils import stats, logger

STUDY = os.environ.get('STUDY')
INPUT_SRC = os.environ.get('INPUT_SRC')
MASIC_PARAM_FILE = os.path.join(INPUT_SRC, "parameters", os.environ.get('MASIC_PARAM_FILE'))
MSGFPLUS_PARAM_FILE= os.path.join(INPUT_SRC, "parameters", os.environ.get('MSGFPLUS_PARAM_FILE'))
MSGFPLUS_MODEF_PARAM_FILE= os.path.join(INPUT_SRC, "parameters", os.environ.get('MSGFPLUS_MODEF_PARAM_FILE'))
MASS_CORRECTION_PARAM= os.path.join(INPUT_SRC, "parameters", os.environ.get('MASS_CORRECTION_PARAM'))

class ProcessingTools:
    '''
    1. .raw -> (msconvert) -> .mzML -> (MSGFPlus) -> .mzid -> (MzidToTsvConverter) -> .tsv -> (TsvToSynConverter) -> _syn.txt
    2. .raw -> (masic) -> _SICstats.txt
    '''

    @stats
    def execute_masic(self, raw_file, raw_basename, final_out, logs):
        '''
         [/I:InputFilePath
         [/O:OutputDirectoryPath] nmdc_jobs/SIC/
         [/P:ParamFilePath]
         [/SF:StatusFileName]
         [/L:[LogFilePath]]
        :param raw_file:
        :param final_out:
        :param logs:
        :return:
        '''

        # MASIC | INPUT:(Thermo .Raw file) | OUTPUT:(_SICStats.txt)
        masic_params = MASIC_PARAM_FILE
        sic = os.path.join(final_out, "SIC")
        if not os.path.exists(sic):
            os.makedirs(sic)
        output_file = os.path.join(sic, raw_basename + "_ScanStats.txt")
        if not os.path.isfile(output_file):
            if os.system("mono /app/masic/MASIC_Console.exe \
                            /I:{} \
                            /O:{} \
                            /P:{} \
                            /SF:{} \
                            /L:{} | tee -a {}".format(raw_file,
                                                      sic,
                                                      masic_params,
                                                      os.path.join(logs, "MasicStatus.xml"),
                                                      os.path.join(logs, "0_masic.commandlog"),
                                                      os.path.join(logs, "0_masic.log")
                                                      )) != 0:
                raise
                # print("Finished running masic")
            # except:
            #     logger.error('MASIC failed for {}-dataset: || RAW: \n {}'.format(STUDY, dataset_id, raw_basename, exc_info=e))
            # except Exception as masic_failed:
            #     raise #masic_error() from masic_failed
        else:
            logger.info("Already exists :: SIC files @:{}".format(sic))

    @stats
    def execute_msconvert(self, raw_file, raw_basename, output, logs):
        '''
        1. MSconvert in pwiz | INPUT:(Thermo .raw file) | OUTPUT:(.mzML file)

        -z [ --zlib ] : use zlib compression for binary data
        --filter arg : add a spectrum list filter
                       peakPicking [<PickerType> [snr=<minimum signal-to-noise ratio>] [peakSpace=<minimum peak spacing>] [msLevel=<ms_levels>]]
                                 : This filter performs centroiding on spectra with the selected <ms_levels>, expressed as an int_set. The value for <PickerType> must be "cwt" or "vendor": when <PickerType> = "vendor", vendor (Windows DLL) code is used if available.
        -o [ --outdir ] arg (=.) : set output directory ('-' for stdout) [.]
        -v [ --verbose ] : display detailed progress information

        --mzML : write mzML format [default]

        :param raw_file:
        :param output: msgfplus_input/
        :param logs:
        :return:
        '''
        # TODO: Only run if mzML doesn't exist
        # if not os.path.exists(output):
        output_file = os.path.join(output, raw_basename + ".mzML")
        if not os.path.isfile(output_file):
            if os.system("wine msconvert \
                    {} \
                    --zlib \
                    --filter 'peakPicking true 2-' \
                    -o {}\
                    --verbose | tee -a {}".format(raw_file,
                                                  output,
                                                  os.path.join(logs, "1_MSconvert.log")
                                                  )) != 0:
                raise
            #     print("Finished running msconvert")
            # except Exception as msconvert_failed:
            #     raise #msconvert_failed()
        else:
            logger.info("Already exists :: mzML file @:{}".format(output))

    @stats
    def execute_msgfplus(self, raw_basename, input, output, fasta, logs):
        '''
        2. MSGFPlus | INPUT:(.mzML file and .fasta file) | OUTPUT:(	.mzid file)
         JVM runs with fixed available memory. Once this memory is exceeded you will receive "java.lang.OutOfMemoryError".
         you can set: https://www.baeldung.com/jvm-parameters

        [-s SpectrumFile] (*.mzML, *.mzXML, *.mgf, *.ms2, *.pkl or *_dta.txt)
        [-o OutputFile]  (*.mzid)
        [-d DatabaseFile] (*.fasta or *.fa or *.faa)
        [-thread NumThreads] (Number of concurrent threads to be executed, Default: Number of available cores)
        [-conf ConfigurationFile]
        [-verbose 0/1] (0: Report total progress only (Default), 1: Report total and per-thread progress/status)
        :param raw_basename:
        :param input:  msgfplus_input/
        :param output: msgfplus_output/
        :param fasta:
        :param logs:
        :return:
        '''
        msgfplus_params = MSGFPLUS_PARAM_FILE
        # FIXME: revCAT location! since, multiple datasets use same FASTA, so they can't generate file in FASTA folder!
        # print('???',input)
        # print('???',raw_basename)
        # if not os.path.exists(output):
        output_file = os.path.join(output, raw_basename + ".mzid")
        if not os.path.isfile(output_file):
            if os.system("java -Xmx32G -jar /app/msgf/MSGFPlus.jar \
                            -s {} \
                            -o {} \
                            -d {} \
                            -thread 16 \
                            -conf {} \
                            -verbose 1 | tee -a {}".format(os.path.join(input, raw_basename + ".mzML"),
                                                           output_file,
                                                           fasta,
                                                           msgfplus_params,
                                                           os.path.join(logs, "2_MSGFPlus.log")
                                                           )) != 0:
                raise
            #     print("Finished running msgfplus")
            # except Exception as msgfplus_failed:
            #     raise #msgfplus_failed()
        else:
            logger.info("Already exists :: mzid file @:{}".format(output))

    @stats
    def execute_mzid2tsv(self, raw_basename, folder, logs):
        '''
        3. MzidToTsvConverter | INPUT:(.mzid file)	| OUTPUT:(.tsv file)
        -mzid:path : Path to the .mzid or .mzid.gz file
        -tsv:path : Path to tsv file to be writte
        -unroll : Signifies that results should be unrolled, giving one line per unique peptide/protein combination in each spectrum identification
        -showDecoy : - Signifies that decoy results should be included in the output .tsv file.
                     - Decoy results have protein names that start with XXX_

        :param raw_basename:
        :param folder: msgfplus_output/
        :param logs:
        :return:
        '''
        tsv_file = os.path.join(folder, raw_basename + ".tsv")
        if not os.path.isfile(tsv_file):
            if os.system("mono /app/mzid2tsv/net462/MzidToTsvConverter.exe \
                            -mzid:{} \
                            -tsv:{} \
                            -unroll \
                            -showDecoy | tee -a {}".format(os.path.join(folder, raw_basename + ".mzid"),
                                                           tsv_file,
                                                           os.path.join(logs, "3_MzidToTsvConverter.log")
                                                           )) != 0:
                raise
            #     print("Finished running MzidToTsvConverter")
            # except Exception as mzid2tsv_failed:
            #     raise #mzid2tsv_failed()
        else:
            logger.info("Already exists :: TSV file @:{}".format(folder))

    @stats
    def execute_tsv2syn(self, raw_basename, input, final_out, revCatfasta, logs):
        '''
        4. TsvToSynConverter| INPUT:(.tsv file) | OUTPUT:(_syn.txt file)

        -I:InputFilePath : MSGF+ results file (_msgfplus.tsv or _msgfdb.tsv or .tsv)
        -O:OutputDirectoryPath
        -M:ModificationDefinitionFilePath
`       -T:MassCorrectionTagsFilePath
        -N:SearchToolParameterFilePath
        -SynPvalue:
        -SynProb:
        -L:LogFilePath
        -ProteinMods to indicate that the _ProteinMods.txt file should be created.
        -F:FastaFilePath]

        :param raw_basename:
        :param input:       msgfplus_output/
        :param final_out:   nmdc_jobs/SYNOPSIS/ [SEQUEST Synopsis/First Hits files]
        :param revCatfasta: Generated by MSGF+
        :param logs:
        :return:
        '''
        msgfplus_params = MASIC_PARAM_FILE
        msgfplus_ModDef_params = MSGFPLUS_MODEF_PARAM_FILE
        mc_params = MASS_CORRECTION_PARAM
        synopsis = os.path.join(final_out, "SYNOPSIS")
        if not os.path.exists(synopsis):
            os.makedirs(synopsis)

        output_file = os.path.join(synopsis, raw_basename + "_syn.txt")
        if not os.path.isfile(output_file):
            if os.system("mono /app/phrp/PeptideHitResultsProcRunner.exe \
                            -I:{} \
                            -O:{} \
                            -M:{} \
                            -T:{} \
                            -N:{} \
                            -SynPvalue:0.2 -SynProb:0.05 \
                            -L:{} \
                            -ProteinMods \
                            -F:{} \
                            | tee -a {}".format(os.path.join(input, raw_basename + ".tsv"),
                                                synopsis,
                                                msgfplus_ModDef_params, mc_params, msgfplus_params,
                                                os.path.join(logs, "4_TsvToSynConverter.commandlog"),
                                                revCatfasta,
                                                os.path.join(logs, "4_TsvToSynConverter.log")
                                                )) != 0:
                raise
            #     print("Finished running TsvToSynConverter")
            # except Exception as tsv2syn_failed:
            #     raise #tsv2syn_failed()
        else:
            logger.info("Already exists :: SYN files @:{}".format(synopsis))

    def convert_fasta_to_txt(self, fasta):
        converted_fasta = os.path.splitext(fasta)[0] + ".txt"
        if not os.path.isfile(converted_fasta):
            if os.system("wine /app/ProteinDigestionSimulator/ProteinDigestionSimulator.exe \
                                -I:{} \
                                -O:{} \
                                -F    \
                                | tee -a {}".format(fasta,
                                                    os.path.split(fasta)[0],
                                                    os.path.join(os.path.split(fasta)[0],
                                                                 "ProteinDigestionSimulator.log"))) != 0:
                raise
        else:
            logger.info("Already exists :: fasta converted files @:{}".format(converted_fasta))

    def get_fasta_loc(self, file_map, dataset_id):
        '''

        :param file_map:
        :param dataset_id:
        :return: path to respective fasta file.
        '''
        # fasta_loc = os.path.join(INPUT_SRC, "fastas/{}".format(STUDY))
        fasta_loc = os.path.join(INPUT_SRC, "fastas/Stegen_IMG_annotations/1781_100336/")
        fptr = open(file_map, "r")
        emsl_to_jgi = json.loads(fptr.read())

        try:
            fasta = emsl_to_jgi[dataset_id][0]
            fasta_path = os.path.join(fasta_loc, *fasta.split("/")[-1:])
            print(">>", fasta_path)
            self.convert_fasta_to_txt(fasta_path)
            return fasta_path
        except KeyError:
            print("Can't find dataset_id:{} on emsl_to_jgi.json file!, so moving on!".format(dataset_id))

    @stats
    def process_datasets(self, data_path):
        '''
        Run dataset through packages.
        make .raw and .fasta files available for them.
        :param data_path:
        :return:
        '''
        file_map = os.path.join(data_path, 'emsl_to_jgi_stegen01152021.json')

        if os.path.isfile(file_map):
            for path, subdirs, files in os.walk(data_path):
                for file in files:
                    if fnmatch.fnmatch(file, '*.raw'):
                        raw_file = os.path.join(path, file)

                        dataset_id = path.split("/")[-1]
                        dataset_faa = self.get_fasta_loc(file_map, dataset_id)
                        if dataset_faa is None:
                            # Can't find an entry in .json.!
                            logger.info(
                                "|---------------FASTA is not available for {}:{} in emsl_to_jgi.json file!---------------".format(
                                    STUDY, dataset_id))
                            continue
                        print('??dataset_faa', dataset_faa)

                        revCat_faa = dataset_faa.replace("faa", "revCat.fasta")
                        print('??revCat_faa', revCat_faa)

                        raw_basename = os.path.basename(os.path.splitext(raw_file)[0])

                        print('``', raw_basename)
                        print('``', raw_file)
                        logs = os.path.join(path, "logs")
                        if not os.path.exists(logs):
                            os.makedirs(logs)

                        ip = os.path.join(path, "msgfplus_input")
                        out = os.path.join(path, "msgfplus_output")
                        nmdc_out = os.path.join(path, "nmdc_jobs")

                        print("processing--------------------------------{}:{}--------------------------------".format(
                            STUDY, dataset_id))
                        # print("Start SELECTED ION CHROMATOGRAMS FILE generation")
                        try:
                            self.execute_masic(raw_file, raw_basename, nmdc_out, logs)
                        except Exception as e:
                            logger.error(
                                'MASIC failed for {}-dataset: || RAW: \n {}'.format(STUDY, dataset_id, raw_basename,
                                                                                    exc_info=e))
                        # print("END SELECTED ION CHROMATOGRAMS FILE generation")
                        # print('*' * 30)
                        # print("Start SYNOPSIS/FIRST-HITS FILE generation")

                        try:
                            self.execute_msconvert(raw_file, raw_basename, ip, logs)
                            try:
                                self.execute_msgfplus(raw_basename, ip, out, dataset_faa, logs)
                                try:
                                    self.execute_mzid2tsv(raw_basename, out, logs)
                                    try:
                                        self.execute_tsv2syn(raw_basename, out, nmdc_out, revCat_faa, logs)
                                        logger.info(
                                            "|--------------------------------{}:{}--------------------------------".format(
                                                STUDY, dataset_id))
                                    except Exception as e:
                                        logger.error(
                                            'TsV2SyN failed for {}-dataset: || RAW: \n {}'.format(STUDY, dataset_id,
                                                                                                  raw_basename,
                                                                                                  exc_info=e))
                                except Exception as e:
                                    logger.error(
                                        'MzID2TsV failed for {}-dataset: || RAW: \n {}'.format(STUDY, dataset_id,
                                                                                               raw_basename,
                                                                                               exc_info=e))
                            except Exception as e:
                                logger.error('MSGFPLUS failed for {}-dataset: || RAW: \n {}'.format(STUDY, dataset_id,
                                                                                                    raw_basename,
                                                                                                    exc_info=e))
                        except Exception as e:
                            logger.error(
                                'MSCONVERT failed for {-/dataset: || RAW: \n {}'.format(STUDY, dataset_id, raw_basename,
                                                                                        exc_info=e))
                        print("End SYNOPSIS/FIRST-HITS FILE generation")
        else:
            print("Can't run without emsl_to_jgi.json? ")


if __name__ == '__main__':

    data = os.path.join(INPUT_SRC, "data/set_of_Dataset_IDs/{}".format(STUDY))

    # TODO shift iteration on datasets here
    tools = ProcessingTools()
    tools.process_datasets(data)


