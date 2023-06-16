import os
import json
import csv
import sys
import hashlib

from urllib.parse import urljoin
from datetime import datetime
from linkml_runtime.dumpers import json_dumper
from nmdc_schema.nmdc import Database, MetaproteomicsAnalysisActivity, DataObject, ProteinQuantification, PeptideQuantification
from nmdc_schema.nmdc_data import get_nmdc_jsonschema_string
from nmdc_id_source import NmdcIdSource


class GenMetadata:
    """
    Generate metadata for the pipeline!
    """
    def __init__(self, analysis_type: str, execution_resource: str, git_url: str, results_url: str,
                  id_source: NmdcIdSource, activity_id: str):

        self.execution_resource = execution_resource
        self.git_url = git_url
        self.type = analysis_type
        self.results_url = results_url
        self.id_source = id_source
        self.activity_id = activity_id

        self.dataset_id = None
        self.genome_directory = None
        self.annotation_file_name = None
        self.resultant_file = None
        self.fasta_file = None
        self.contaminant_file = None
        self.peptide_report = None
        self.protein_report = None
        self.qc_metric_report = None
        self.started_at_time = None
        self.ended_at_time = None

    def get_md5(self, file):
        """
        generates MD5 checksum of a file.

        :param file: filename with absolute path.
        :return: chechsum or empty string.
        """
        if not file in (None, ""):
            md5 = hashlib.md5(open(file, "rb").read()).hexdigest()
            return md5
        else:
            return ""

    def get_ProteinQuantification(self):
        protein_arr = []

        with open(self.protein_report, encoding="utf-8-sig") as tsvfile:
            tsvreader = csv.reader(tsvfile, delimiter="\t")
            header = next(tsvreader)

            for line in tsvreader:
                protein_quantification = ProteinQuantification(
                    peptide_sequence_count = line[1],
                    best_protein = line[2].replace(" ", ""),
                    all_proteins = [protein.replace(" ", "") for protein in line[9].split(",")],
                    protein_spectral_count = line[12],
                    protein_sum_masic_abundance = line[13],
                )

                protein_arr.append(protein_quantification)

    def get_PeptideQuantification(self):
        peptide_quantifications_arr = []

        with open(self.peptide_report, encoding="utf-8-sig") as tsvfile:
            tsvreader = csv.reader(tsvfile, delimiter="\t")
            header = next(tsvreader)
            for line in tsvreader:
                peptide_quantification = PeptideQuantification(
                peptide_sequence = line[1],
                best_protein = line[2].replace(" ", ""),
                all_proteins = [protein.replace(" ", "") for protein in line[9].split(",")],
                min_q_value = line[11],
                peptide_spectral_count = line[12],
                peptide_sum_masic_abundance = int(float(line[13]))
                )

                peptide_quantifications_arr.append(peptide_quantification)
        
        return peptide_quantifications_arr

    def get_file_data_object(self, file_path, file_name, description, activity_id):
        """
        - Makes entry in _emsl_analysis_data_objects.json
        - Contains information about Pipeline's analysis results E.g.
            1. resultant.tsv
            2. data_out_table.tsv
        """

        file_id = self.id_source.get_ids("nmdc:DataObject", 1)[0]
        
        data_object = DataObject(
            id=file_id,
            name=file_name,
            description=description,
            file_size_bytes=os.stat(file_path).st_size,
            md5_checksum=hashlib.md5(open(file_path,'rb').read()).hexdigest(),
            was_generated_by=activity_id,
        )

        if "MSGFjobs_MASIC_resultant" in file_name:
            data_object.data_object_type = "Unfiltered Metaproteomics Results"
        if "Peptide_Report" in file_name:
            data_object.data_object_type = "Peptide Report"
        if "Protein_Report" in file_name:
            data_object.data_object_type = "Protein Report"
        if "QC_metrics" in file_name:
            data_object.data_object_type = "Metaproteomics Workflow Statistics"

        data_object.url = urljoin(self.results_url, file_name)
        data_object.was_generated_by = self.activity_id

        return data_object


    def get_data_objects(self):

        data_objects = []

        data_objects.append(
            self.get_file_data_object(
                self.resultant_file,
                os.path.basename(self.resultant_file),
                "Aggregation of analysis tools {MSGFplus, MASIC} results", 
                self.activity_id
            )
        )

        data_objects.append(
            self.get_file_data_object(
                self.peptide_report,
                os.path.basename(self.peptide_report),
                "Aggregated peptides sequences from MSGF+ search results filtered to ~5% FDR", 
                self.activity_id
            )
        )

        data_objects.append(
            self.get_file_data_object(
                self.protein_report,
                os.path.basename(self.protein_report),
                "Aggregated protein lists from MSGF+ search results filtered to ~5% FDR", 
                self.activity_id
            )
        )

        data_objects.append(
            self.get_file_data_object(
                self.qc_metric_report,
                os.path.basename(self.qc_metric_report),
                "Overall statistics from MSGF+ search results filtered to ~5% FDR", 
                self.activity_id
            )
        )

        return data_objects

    def get_has_input(self) -> list:
        has_input = []

        # add .RAW
        raw_file_name = self.dataset_id
        has_input.append(raw_file_name)
        
        fasta_checksum = self.get_md5(self.fasta_file)
        if fasta_checksum:
            # add to metadata.
            has_input.append("nmdc:" + fasta_checksum)
        else:
            print("Found HASH empty for {}".format(self.fasta_file))

        # add EMSL contaminants
        contaminant_checksum = self.get_md5(self.contaminant_file)
        if contaminant_checksum:
            has_input.append("nmdc:" + contaminant_checksum)
        else:
            print("Found HASH empty for {}".format(self.fasta_file))
        
        return has_input

    def get_metaproteomics_analysis_activity(self):
        has_input_arr = self.get_has_input()
        has_output_arr = False

        mp_analysis_activity_obj = MetaproteomicsAnalysisActivity(
            id=self.activity_id,
            execution_resource=self.execution_resource,
            git_url=self.git_url,
            name=":".join(["Metaproteome", self.dataset_id, self.genome_directory]),
            was_informed_by=":".join(["emsl", self.dataset_id]),
            type=self.type,
            has_output=has_output_arr,
            has_input=has_input_arr,
            started_at_time="2002-09-24-06:00",
            ended_at_time="2002-09-24-06:00"
            )

        mp_analysis_activity_obj.has_peptide_quantifications = self.get_PeptideQuantification()

        return mp_analysis_activity_obj
    
    def set_keys(
        self,
        dataset_id,
        genome_directory,
        annotation_file_name,
        resultant_file,
        fasta_file,
        contaminant_file,
        peptide_report,
        protein_report,
        qc_metric_report,
        started_at_time,
        ended_at_time,
    ):
        self.dataset_id = dataset_id
        self.genome_directory = genome_directory
        self.annotation_file_name = annotation_file_name
        self.resultant_file = resultant_file
        self.fasta_file = fasta_file
        self.contaminant_file = contaminant_file
        self.peptide_report = peptide_report
        self.protein_report = protein_report
        self.qc_metric_report = qc_metric_report
        self.started_at_time = started_at_time
        self.ended_at_time = ended_at_time

    @staticmethod
    def create(analysis_type: str, execution_resource: str, git_url: str, results_url: str,
                  id_source: NmdcIdSource) -> 'GenMetadata':
        '''
        Create a GenMetadata object with newly minted activity ID with OOP in mind.
        '''
        activity_id = id_source.get_ids(analysis_type, 1)[0]

        return GenMetadata(analysis_type, execution_resource, git_url, results_url, id_source, activity_id)


if __name__ == "__main__":
    mapper_file = sys.argv[1]
    study = sys.argv[2]
    execution_resource = sys.argv[3]
    git_url = sys.argv[4]
    results_url = sys.argv[5]

    if results_url[-1] != '/':
        results_url = results_url + '/'

    # An array of DataObject
    data_objects_arr = []
    # An array of ""
    metaproteomics_analysis_activity_arr = []

    # Minting source
    # TODO grab client_id and client_secret from env
    id_source = NmdcIdSource("", "")

    # 1. Make collection and populate them.
    with open(mapper_file, "r+") as json_file:
        mapper_json = json.load(json_file)

        for mapping in mapper_json:

            meta_file = GenMetadata.create(
                "nmdc:MetaproteomicsAnalysisActivity",
                execution_resource=execution_resource,
                git_url=git_url,
                results_url=results_url,
                id_source=id_source
                )

            meta_file.set_keys(
                mapping["dataset_id"],
                mapping["genome_directory"],
                os.path.basename(mapping["txt_faa_file"]),
                mapping["resultant_file"],
                mapping["faa_file"],
                mapping["contaminate_file"],
                mapping["peptide_report_file"],
                mapping["protein_report_file"],
                mapping["qc_metric_report_file"],
                mapping["start_time"],
                mapping["end_time"],
            )

            metaproteomics_analysis_activity = meta_file.get_metaproteomics_analysis_activity()
            data_objects = meta_file.get_data_objects()

            metaproteomics_analysis_activity_arr.append(metaproteomics_analysis_activity)
            data_objects_arr.extend(data_objects)
            

        # 2. dump collections in json files
        schema = get_nmdc_jsonschema_string()
        activity_file = f"{study}_MetaProteomicAnalysis_activity.json"
        data_obj_file = f"{study}_analysis_data_objects.json"

        activity_db = Database()
        activity_db.metaproteomics_analysis_activity_set = metaproteomics_analysis_activity_arr

        data_obj_db = Database()
        data_obj_db.data_object_set = data_objects_arr

        json_dumper.dump(activity_db, contexts=schema, inject_type=False, to_file=activity_file)
        json_dumper.dump(data_obj_db, contexts=schema, inject_type=False, to_file=data_obj_file)