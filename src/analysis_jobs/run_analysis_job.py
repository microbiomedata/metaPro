import os
import json
from datetime import datetime
import time
import tracemalloc
import copy
from utility.utils import logger

from src.analysis_jobs.merge_jobs.DatasetsMerger import DatasetsMerger

LOGGED_ANALYSIS_JOB = {}
ANALYSIS_JOBS_OBJECT = {}


def stats(method):
    def collected(*args, **kw):
        ts = time.time()
        tracemalloc.start()
        result = method(*args, **kw)
        current, peak = tracemalloc.get_traced_memory()
        te = time.time()
        if method.__name__ == "analysis_job":
            if "total_runtime" not in ANALYSIS_JOBS_OBJECT:
                ANALYSIS_JOBS_OBJECT["total_runtime"] = {
                    "runtime": time.strftime("%H:%M:%S", time.gmtime((te - ts))),
                    "current_m/m(MB)": current / 10 ** 6,
                    "peak_m/m(MB)": peak / 10 ** 6,
                }
        else:
            if method.__name__ not in ANALYSIS_JOBS_OBJECT:
                ANALYSIS_JOBS_OBJECT[method.__name__] = {
                    "runtime": time.strftime("%H:%M:%S", time.gmtime((te - ts))),
                    "current_m/m(MB)": current / 10 ** 6,
                    "peak_m/m(MB)": peak / 10 ** 6,
                }
        tracemalloc.stop()
        return result

    return collected


class ProcessingTools:
    """
    1. .raw -> (masic) -> _SICstats.txt
    2. .raw -> (msconvert) -> .mzML -> (MSGFPlus) -> .mzid -> (MzidToTsvConverter) -> .tsv -> (TsvToSynConverter) -> _syn.txt
    """

    def __init__(self, result_loc=None, mappings=None, parameters=None):
        self.result_loc = result_loc
        self.save_job_results = None
        self.mappings = mappings
        self.logs_collected_at = None
        self.dataset_name = None
        self.contaminant_file_loc = None

        self.MASIC_PARAM_FILE = parameters["MASIC_PARAM_FILE"]
        self.MSGFPLUS_PARAM_FILE = parameters["MSGFPLUS_PARAM_FILE"]
        self.MSGFPLUS_MODEF_PARAM_FILE = parameters["MSGFPLUS_MODEF_PARAM_FILE"]
        self.MASS_CORRECTION_PARAM_FILE = parameters["MASS_CORRECTION_PARAM_FILE"]

    def register_job_in_emsl_to_jgi(
        self, dataset_id, genome_directory, key, value, emsl_to_jgi_copy
    ):

        locations = emsl_to_jgi_copy[dataset_id]["genome_directory"][genome_directory]
        if key not in locations:
            locations[key] = value
        else:
            print(f"{key} already present in {genome_directory}")

    @stats
    def execute_masic(self, raw_file):
        """
         [/I:InputFilePath
         [/O:OutputDirectoryPath] nmdc_jobs/SIC/
         [/P:ParamFilePath]
         [/SF:StatusFileName]
         [/L:[LogFilePath]]
        :param raw_file:
        :return:
        """
        save_at = os.path.join(self.save_job_results, "nmdc_jobs", "SIC")
        if not os.path.exists(save_at):
            os.makedirs(save_at)

        output_file = os.path.join(save_at, self.dataset_name + "_ScanStats.txt")
        if not os.path.isfile(output_file):
            if (
                os.system(
                    f"mono /app/masic/MASIC_Console.exe \
                            /I:{raw_file} \
                            /O:{save_at} \
                            /P:{self.MASIC_PARAM_FILE} \
                            /SF:{ os.path.join( self.log_collected_at, 'MasicStatus.xml') } \
                            /L:{os.path.join( os.path.abspath(self.save_job_results), 'analysis_jobs_logs',  '0_masic.commandlog') }  \
                            | tee -a { os.path.join( self.log_collected_at,'0_masic.log') }"
                )
                != 0
            ):
                raise
        else:
            print(f"Already exists :: {output_file}")

    @stats
    def execute_msconvert(self, raw_file):
        """
        1. MSconvert in pwiz | INPUT:(Thermo .raw file) | OUTPUT:(.mzML file)

        -z [ --zlib ] : use zlib compression for binary data
        --filter arg : add a spectrum list filter
                       peakPicking [<PickerType> [snr=<minimum signal-to-noise ratio>] [peakSpace=<minimum peak spacing>] [msLevel=<ms_levels>]]
                                 : This filter performs centroiding on spectra with the selected <ms_levels>, expressed as an int_set. The value for <PickerType> must be "cwt" or "vendor": when <PickerType> = "vendor", vendor (Windows DLL) code is used if available.
        -o [ --outdir ] arg (=.) : set output directory ('-' for stdout) [.]
        -v [ --verbose ] : display detailed progress information

        --mzML : write mzML format [default]
        :param raw_file:
        :return:
        """
        save_at = os.path.join(self.save_job_results, "msgfplus_input")
        if not os.path.exists(save_at):
            os.makedirs(save_at)

        output_file = os.path.join(save_at, self.dataset_name + ".mzML")
        if not os.path.isfile(output_file):
            if (
                os.system(
                    f"wine msconvert \
                                {raw_file} \
                                --zlib \
                                --filter 'peakPicking true 2-' \
                                -o {save_at}\
                                --verbose | tee -a {os.path.join( self.log_collected_at, '1_MSconvert.log')}"
                )
                != 0
            ):
                raise
        else:
            print(f"Already exists :: {output_file}")
        return output_file

    @stats
    def execute_msgfplus(self, mzml_file, fasta_file):
        """
        2. MSGFPlus | INPUT:(.mzML file and .fasta file) | OUTPUT:(	.mzid file)
         JVM runs with fixed available memory. Once this memory is exceeded you will receive "java.lang.OutOfMemoryError".
         you can set: https://www.baeldung.com/jvm-parameters

        [-s SpectrumFile] (*.mzML, *.mzXML, *.mgf, *.ms2, *.pkl or *_dta.txt)
        [-o OutputFile]  (*.mzid)
        [-d DatabaseFile] (*.fasta or *.fa or *.faa)
        [-thread NumThreads] (Number of concurrent threads to be executed, Default: Number of available cores)
        [-conf ConfigurationFile]
        [-verbose 0/1] (0: Report total progress only (Default), 1: Report total and per-thread progress/status)

        :param mzml_file:
        :param fasta_file: contaminated fasta file.
        :return:
        """
        save_at = os.path.join(self.save_job_results, "msgfplus_output")
        if not os.path.exists(save_at):
            os.makedirs(save_at)

        output_file = os.path.join(save_at, self.dataset_name + ".mzid")
        if not os.path.isfile(output_file):
            if (
                os.system(
                    f"java -Xmx32G -jar /app/msgf/MSGFPlus.jar \
                            -s {mzml_file} \
                            -o {output_file} \
                            -d {fasta_file} \
                            -thread 16 \
                            -conf {self.MSGFPLUS_PARAM_FILE} \
                            -verbose 1 | tee -a {os.path.join( self.log_collected_at, '2_MSGFPlus.log')}"
                )
                != 0
            ):
                raise

        else:
            print(f"Already exists :: {output_file}")

        # return the path to new fasta loc to run the tsvtosyn!
        faa_name = os.path.splitext((os.path.basename(fasta_file)))[0]
        revcat_fasta = os.path.join(
            save_at, "fasta_residuals", f"{faa_name}.revCat.fasta"
        )
        return output_file, revcat_fasta

    @stats
    def execute_mzid2tsv(self, mzid_file):
        """
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
        """
        save_at = os.path.join(self.save_job_results, "msgfplus_output")
        if not os.path.exists(save_at):
            os.makedirs(save_at)

        output_file = os.path.join(save_at, self.dataset_name + ".tsv")
        if not os.path.isfile(output_file):
            if (
                os.system(
                    f"mono /app/mzid2tsv/net462/MzidToTsvConverter.exe \
                            -mzid:{mzid_file} \
                            -tsv:{output_file} \
                            -unroll \
                            -showDecoy | tee -a {os.path.join( self.log_collected_at, '3_MzidToTsvConverter.log')}"
                )
                != 0
            ):
                raise
        else:
            print(f"Already exists :: {output_file}")
        return output_file

    @stats
    def execute_tsv2syn(self, tsv_file, revcat_faa_file):
        """
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
                :param tsv_file:
                :param revcat_faa_file:
                :return:
        """
        save_at = os.path.join(self.save_job_results, "nmdc_jobs", "SYNOPSIS")
        if not os.path.exists(save_at):
            os.makedirs(save_at)

        output_file = os.path.join(save_at, self.dataset_name + "_syn.txt")
        if not os.path.isfile(output_file):
            if (
                os.system(
                    f"mono /app/phrp/PeptideHitResultsProcRunner.exe \
                            -I:{tsv_file} \
                            -O:{save_at} \
                            -M:{self.MSGFPLUS_MODEF_PARAM_FILE} \
                            -T:{self.MASS_CORRECTION_PARAM_FILE} \
                            -N:{self.MSGFPLUS_PARAM_FILE} \
                            -SynPvalue:0.2 -SynProb:0.05 \
                            -L:{ os.path.join( self.log_collected_at, '4_TsvToSynConverter.commandlog')} \
                            -ProteinMods \
                            -F:{revcat_faa_file} \
                            | tee -a {os.path.join( self.log_collected_at, '4_TsvToSynConverter.log')}"
                )
                != 0
            ):
                raise
        else:
            print(f"Already exists :: {output_file}")

    @stats
    def execute_protein_digestion_simulator(self, fasta_file):
        """
        Runs on fasta_file+contaiminant_file
        :param fasta_file:
        :return:
        """
        save_at = os.path.join(self.save_job_results, "protein_digestion")
        if not os.path.exists(save_at):
            os.makedirs(save_at)

        output_file = os.path.join(
            save_at, os.path.splitext((os.path.basename(fasta_file)))[0] + ".txt"
        )
        if not os.path.isfile(output_file):
            if (
                os.system(
                    f"wine /app/ProteinDigestionSimulator/ProteinDigestionSimulator.exe \
                                -I:{fasta_file} \
                                -O:{save_at} \
                                -F    \
                                | tee -a { os.path.join( self.log_collected_at, 'ProteinDigestionSimulator.log') }"
                )
                != 0
            ):

                raise
        else:
            print(f"Already exists :: {output_file}")

        return output_file

    @stats
    def analysis_job(self, dataset_id, faa_file_loc, raw_file_loc):
        """

        :param dataset_id:
        :param faa_file_loc: contaminated fasta file.
        :param raw_file_loc:
        :return:
        """
        print(f"Processing start {dataset_id} : {os.path.basename(faa_file_loc)}")
        # 1..raw -> (masic) -> _SICstats.txt
        try:
            self.execute_masic(raw_file_loc)
        except Exception as e:
            print(f"MASIC_Console failed for {dataset_id}: \n {e}")
        # 2. .raw -> (msconvert) -> .mzML -> (MSGFPlus) -> .mzid -> (MzidToTsvConverter) -> .tsv -> (TsvToSynConverter) -> _syn.txt

        try:
            mzml_file = self.execute_msconvert(raw_file_loc)
            try:
                mzid_file, revcat_faa_file = self.execute_msgfplus(
                    mzml_file, faa_file_loc
                )
                try:
                    tsv_file = self.execute_mzid2tsv(mzid_file)
                    try:
                        self.execute_tsv2syn(tsv_file, revcat_faa_file)
                        print(
                            f"Processing ended {dataset_id} : {os.path.basename(faa_file_loc)}"
                        )

                    except Exception as e:
                        logger.error(
                            f"PeptideHitResultsProcRunner failed for {dataset_id}",
                            exc_info=True,
                        )
                except Exception as e:
                    logger.error(
                        f"MzidToTsvConverter failed for {dataset_id}", exc_info=True
                    )
            except Exception as e:
                logger.error(f"MSGFPlus failed for {dataset_id}", exc_info=True)
        except Exception as e:
            logger.error(f"msconvert failed for {dataset_id}", exc_info=True)
        pass

    def run_n_log_job(
        self, dataset_id, genome_directory, faa_file_loc, raw_file_loc, emsl_to_jgi_copy
    ):

        job_name = f"{dataset_id}:{genome_directory}"
        job_object = {"job_name": job_name}
        self.register_job_in_emsl_to_jgi(
            dataset_id, genome_directory, "job_name", job_name, emsl_to_jgi_copy
        )
        # log dataset_id
        if dataset_id not in LOGGED_ANALYSIS_JOB:
            LOGGED_ANALYSIS_JOB[dataset_id] = [job_object]
        else:
            LOGGED_ANALYSIS_JOB[dataset_id].append(job_object)

        # log starttime
        started_at_time = datetime.utcnow().strftime("%Y_%m_%d-%I_%M_%S_%p")
        for job in LOGGED_ANALYSIS_JOB[dataset_id]:
            if job["job_name"] == job_name:
                if "started_at_time" not in job:
                    job["started_at_time"] = started_at_time
        self.register_job_in_emsl_to_jgi(
            dataset_id,
            genome_directory,
            "started_at_time",
            started_at_time,
            emsl_to_jgi_copy,
        )

        # run the job
        self.analysis_job(dataset_id, faa_file_loc, raw_file_loc)

        # log endtime
        ended_at_time = datetime.utcnow().strftime("%Y_%m_%d-%I_%M_%S_%p")
        for job in LOGGED_ANALYSIS_JOB[dataset_id]:
            if job["job_name"] == job_name:
                if "ended_at_time" not in job:
                    job["ended_at_time"] = ended_at_time
        self.register_job_in_emsl_to_jgi(
            dataset_id,
            genome_directory,
            "ended_at_time",
            ended_at_time,
            emsl_to_jgi_copy,
        )

        # log analysis_job object
        for job in LOGGED_ANALYSIS_JOB[dataset_id]:
            if job["job_name"] == job_name:
                if "analysis_jobs" not in job:
                    job["analysis_jobs"] = ANALYSIS_JOBS_OBJECT
        self.register_job_in_emsl_to_jgi(
            dataset_id,
            genome_directory,
            "analysis_jobs",
            ANALYSIS_JOBS_OBJECT,
            emsl_to_jgi_copy,
        )
        pass

    def convert_faa2txt(self, dataset_id, faa_file_loc):
        faa_txt = ""
        try:
            faa_txt = self.execute_protein_digestion_simulator(faa_file_loc)
        except Exception as e:
            logger.error(
                f"ProteinDigestionSimulator failed for {dataset_id}", exc_info=True
            )
        return faa_txt

    def merge_analysis_jobs(self, dataset_id, genome_directory):

        save_at = os.path.join(self.save_job_results, "merged_jobs")
        if not os.path.exists(save_at):
            os.makedirs(save_at)

        merge = DatasetsMerger(self.save_job_results)
        resultant_df = merge.merge_all_jobs_in_UserInput()

        output_file = os.path.join(
            save_at, f"{dataset_id}_{genome_directory}_MSGFjobs_MASIC_resultant.tsv"
        )
        if not os.path.isfile(output_file):
            try:
                resultant_df.to_csv(output_file, sep="\t")
            except Exception as e:
                print("Error : " + e)
        return output_file

    def contaminate_fasta(self, files):

        save_at = os.path.join(
            self.save_job_results, "msgfplus_output", "fasta_residuals"
        )
        if not os.path.exists(save_at):
            os.makedirs(save_at)

        contaminated_file = os.path.join(save_at, os.path.basename(files[0]))
        with open(contaminated_file, "w") as outfile:
            for file in files:
                with open(file) as infile:
                    for line in infile:
                        outfile.write(line)
        return contaminated_file

    def process_datasets(self):
        """
        Process raw_file and fasta_file combination from given mappings.
        :param data_path:
        :return:
        """

        with open(self.mappings, "r+") as json_file:
            emsl_to_jgi = json.load(json_file)
            emsl_to_jgi_copy = copy.deepcopy(emsl_to_jgi)

            contaminant_file_loc = emsl_to_jgi["contaminant_file_loc"]
            # run for each dataset
            for dataset_id, values in emsl_to_jgi.items():
                if dataset_id not in [
                    "contaminant_file_loc",
                    "analysis_activity_file_loc",
                    "data_objects_file_loc",
                    "STUDY",
                    "tools_used",
                ]:
                    raw_file_loc = values["raw_file_loc"]
                    self.dataset_name = values["dataset_name"]
                    # dataset search against a fasta file
                    for genome_directory, locations in values[
                        "genome_directory"
                    ].items():
                        # clear object to prepare next job
                        ANALYSIS_JOBS_OBJECT.clear()

                        # create log_dir
                        self.save_job_results = os.path.join(
                            self.result_loc, dataset_id, genome_directory
                        )
                        self.log_collected_at = os.path.join(
                            os.path.abspath(self.save_job_results), "analysis_jobs_logs"
                        )
                        if not os.path.exists(self.log_collected_at):
                            os.makedirs(self.log_collected_at)

                        files = [locations["faa_file_loc"], contaminant_file_loc]
                        contaminated_faa_file_loc = self.contaminate_fasta(files)

                        self.register_job_in_emsl_to_jgi(
                            dataset_id,
                            genome_directory,
                            "contaminated_faa_file_loc",
                            contaminated_faa_file_loc,
                            emsl_to_jgi_copy,
                        )
                        # convert .faa to .txt
                        faa_txt_file = self.convert_faa2txt(
                            dataset_id, contaminated_faa_file_loc
                        )
                        self.register_job_in_emsl_to_jgi(
                            dataset_id,
                            genome_directory,
                            "txt_faa_file_loc",
                            faa_txt_file,
                            emsl_to_jgi_copy,
                        )

                        # log & run job
                        self.run_n_log_job(
                            dataset_id,
                            genome_directory,
                            contaminated_faa_file_loc,
                            raw_file_loc,
                            emsl_to_jgi_copy,
                        )

                        # merge analysis
                        resultant_file = self.merge_analysis_jobs(
                            dataset_id, genome_directory
                        )
                        self.register_job_in_emsl_to_jgi(
                            dataset_id,
                            genome_directory,
                            "resultant_file_loc",
                            resultant_file,
                            emsl_to_jgi_copy,
                        )

            # capture the job metadata object
            logger.info("Jobrun", extra=LOGGED_ANALYSIS_JOB)

            # update emsl_to_jgi.json
            json_file.seek(0)  # move back to BOF.
            json_file.truncate()
            json_file.write(json.dumps(emsl_to_jgi_copy, default=str, indent=4))
        pass

if __name__ == "__main__":

    result_loc = os.path.join("storage", "results", os.environ.get("STUDY"))
    mapper_file = os.path.join(result_loc, "emsl_to_jgi.json")
    param_file_loc = os.path.join("storage", "parameters")
    parameters = {
        "MASIC_PARAM_FILE": os.path.join(
            param_file_loc, os.environ.get("MASIC_PARAM_FILENAME")
        ),
        "MSGFPLUS_PARAM_FILE": os.path.join(
            param_file_loc, os.environ.get("MSGFPLUS_PARAM_FILENAME")
        ),
        "MSGFPLUS_MODEF_PARAM_FILE": (
            os.path.join(
                param_file_loc, os.environ.get("MSGFPLUS_MODEF_PARAM_FILENAME")
            )
        )
        if os.environ.get("MSGFPLUS_MODEF_PARAM_FILENAME") is not None
        else "",
        "MASS_CORRECTION_PARAM_FILE": (
            os.path.join(
                param_file_loc, os.environ.get("MASS_CORRECTION_PARAM_FILENAME")
            )
        )
        if os.environ.get("MSGFPLUS_MODEF_PARAM_FILENAME") is not None
        else "",
    }

    if os.path.isfile(mapper_file):
        tools = ProcessingTools(result_loc, mapper_file, parameters)
        tools.process_datasets()
    else:
        print("Can't run workflow without emsl_to_jgi.json.")
