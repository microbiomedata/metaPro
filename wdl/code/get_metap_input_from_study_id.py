from nmdc_api.schema_api import *
from nmdc_api.utils import download_file
import click
import hashlib
import logging
import csv
from logging import StreamHandler
from logging.handlers import RotatingFileHandler
from pathlib import Path
from pydantic import BaseModel, Field, TypeAdapter
from typing import List, Iterable, Callable, TypeVar, Optional, Tuple, Any
from operator import itemgetter
from dataclasses import dataclass, fields, asdict


def setup_logging(level=logging.INFO, logfile="get_metap_input_from_study_id.log"):
    fmt = "%(asctime)s %(levelname)s %(name)s: %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"

    console = StreamHandler()
    filerot = RotatingFileHandler(logfile, maxBytes=5_000_000, backupCount=3, encoding="utf-8")

    logging.basicConfig(
        level=level,
        format=fmt,
        datefmt=datefmt,
        handlers=[console, filerot],
        force=True,
    )


datafiles_dir = Path("datafiles")
settings_dir = Path("settings")


class DataFilesMapping(BaseModel):
    data_generation_id: str = ""
    dataset_name: str = ""
    raw_file_loc: Path = datafiles_dir
    dataset_id: str = ""
    faa_file_loc: Path = datafiles_dir
    faa_file_id: str = ""
    gff_file_loc: Path = datafiles_dir
    gff_file_id: str = ""
    analysis_id: str = ""


mapping_adapter = TypeAdapter(List[DataFilesMapping])


class MetapInput(BaseModel):
    mapper_list: List[DataFilesMapping] = Field(serialization_alias="metapro.mapper_list", default_factory=list)
    masic_param_file_filepath: Path = Field(serialization_alias="metapro.MASIC_PARAM_FILE_LOC", default=settings_dir)
    msgf_param_file_filepath: Path = Field(serialization_alias="metapro.MSGFPLUS_PARAM_FILE_LOC", default=settings_dir)
    kaiko_param_file_filepath: Path = Field(serialization_alias="metapro.KAIKO_PARAM_FILE_LOC", default=settings_dir)
    contaminant_param_file_filepath: Path = Field(serialization_alias="metapro.CONTAMINANT_FILE_LOC", default=settings_dir)
    qvalue_threshold: str = Field(serialization_alias="metapro.QVALUE_THRESHOLD", default="")
    study: str = Field(serialization_alias="metapro.STUDY", default="")
    execution_resource: str = Field(serialization_alias="metapro.EXECUTION_RESOURCE", default="")
    data_url: str = Field(serialization_alias="metapro.DATA_URL", default="")
    masic_parameter_file_id: str = Field(serialization_alias="metapro.MASIC_PARAM_FILE_ID", default="")
    msgf_parameter_file_id: str = Field(serialization_alias="metapro.MSGFPLUS_PARAM_FILE_ID", default="")
    contaminant_parameter_file_id: str = Field(serialization_alias="metapro.CONTAMINANT_FILE_ID", default="")
    kaiko_param_file_filepath: Path = Field(serialization_alias="metapro.KAIKO_PARAM_FILE_LOC", default=settings_dir)
    metagenome_free: bool = Field(serialization_alias="metapro.METAGENOME_FREE", default=False)


@dataclass 
class DownloadedFile:
    filename: str = ""
    filepath: str = ""
    is_success: bool = False
    do_md5: str = ""
    calculated_md5: str = ""


def colon_to_underscore(s: str) -> str:
    return s.replace(':', '_')


T = TypeVar("T")
def first_or_default(it: Iterable[T],
                     predicate: Callable[[T], bool] | None = None,
                     default: Optional[T] = None) -> Optional[T]:
    if predicate is None:
        return next(iter(it), default)
    return next((x for x in it if predicate(x)), default)


def get_md5(path: Path):
    md5 = hashlib.md5()

    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            md5.update(chunk)

    return md5.hexdigest()


def download_files(colored_do: List[Tuple[str, Any]], settings_dir: Path, datafiles_dir: Path, check_md5: bool) -> List[DownloadedFile]:
    dls: List[DownloadedFile] = []
    for do in colored_do:
        dat = {}
        path = settings_dir if do[0] == "settings" else datafiles_dir

        log.info(f"downloading {do[1]['name']} to {path} from {do[1]['url']}")
        filepath = download_file(do[1]["url"], path)
        is_sucess = filepath is not None
        log.info(f"{'successfully downloaded' if is_sucess else 'failed to download'} {do[1]['name']}")

        if is_sucess:
            if check_md5:
                fp = path / do[1]["name"]
                file_md5 = get_md5(fp)
                is_sucess = do[1]["md5_checksum"] == file_md5
                log.info(f"{'matching' if is_sucess else 'unmatched'} md5 checksum for {do[1]['name']}")
                
                dat["do_md5"] = do[1]["md5_checksum"]
                dat["calculated_md5"] = file_md5
            dat["filepath"] = filepath

        dat["filename"] = do["name"]
        dat["is_success"] = is_sucess
        dls.append(DownloadedFile(**dat))
    
    return dls


setup_logging(logging.DEBUG)
log = logging.getLogger(__name__)

@click.command()
@click.option('--study_id', prompt='Your study id', help='The study id to get biosample IDs from.')
@click.option('--output_dir', prompt='Output directory', help='The absolute output directory to write downloaded files to.')
@click.option('-dl', is_flag=True, default=False, help='Set to download files.')
@click.option('-md5', is_flag=True, default=False, help='Checks md5 of downloaded file if dl is specified.')
def main(study_id, output_dir, dl, md5):

    to_download_list: List[Tuple[str, Any]] = [] 

    # ideally, these would be passed in or specified some other way
    # but so far these are static
    masic_param_id = "nmdc:dobj-11-hfx93f93"
    msgf_param_id = "nmdc:dobj-11-h9637w90"
    contam_id =  "nmdc:dobj-11-sprrem27"
    q_value_threshold = "0.05"
    data_url = "https://nmdcdemo.emsl.pnnl.gov/proteomics/results/"
    execution_resource = "EMSL"

    # set-up paths
    output_path = Path(output_dir)
    settings_path = output_path / "settings"
    datafiles_path = output_path / "datafiles"
    output_path.mkdir(parents=True, exist_ok=True)
    settings_path.mkdir(parents=True, exist_ok=True)
    datafiles_path.mkdir(parents=True, exist_ok=True)

    log.info(f"getting data generation records for {study_id}")

    # get data generation records for study
    set = NMDCCollection.data_generation_set
    filters = [
        filter_field_value_equals("associated_studies", study_id),
    ]

    records = get_records(set, filters)
    if records is None or len(records) == 0:
        log.warning(f"No records found in '{set}' for study id '{study_id}'")
        return
    
    log.info(f"return {len(records)} data generation records for {study_id}")

    grouped_records = {x["has_input"][0]: {} for x in records }

    # add all metagenome data gen records
    for record in records:
        if record["analyte_category"] == "metagenome":
            grouped_records[record["has_input"][0]]["metagenome"] = record

    # walk up processed sample -> material processing to reach parent biosample for association and add 
    # proteome data gen records mapped to biosample instead of processed sample
    for record in records:
        if record["analyte_category"] == "metaproteome":
            set = NMDCCollection.processed_sample_set
            filters = [ filter_field_value_matches_any("id", record["has_input"]) ]
            ret = get_records(set, filters)
            while len(ret) != 0:
                if set == NMDCCollection.material_processing_set:
                    if "has_input" in ret[0] and ret[0]["has_input"][0].startswith("nmdc:bsm"):
                        grouped_records[ret[0]["has_input"][0]]["metaproteome"] = record
                        break
                    set = NMDCCollection.processed_sample_set
                    filters = [ filter_field_value_matches_any("id", ret[0]["has_input"]) ]
                    ret = get_records(set, filters)
                else:
                    set = NMDCCollection.material_processing_set
                    filters = [ filter_field_value_matches_any("has_output", [ ret[0]["id"] ]) ]
                    ret = get_records(set, filters)

    # remove everything that's not paired
    paired_dict = {k: grouped_records[k] for k in grouped_records.keys() if len(grouped_records[k]) == 2}  
    
    # for convenience, build the mapping of (biosample id, dgns id, dgms id)
    dga_map = {
        (k, v["metagenome"]["id"], v["metaproteome"]["id"]): v
        for k,v
        in paired_dict.items()
    }

    # get annotation records from paired dict and download files, build mapping list
    data_files_mapping: List[DataFilesMapping] = []
    to_filter_on = []
    for k, _ in paired_dict.items():
        to_filter_on.append(paired_dict[k]["metagenome"]["id"])

    set = NMDCCollection.workflow_execution_set        
    filters = [
        filter_field_value_equals("type", "nmdc:MetagenomeAnnotation"),
        filter_field_value_matches_any("was_informed_by", to_filter_on)
    ]
    annotation_records = get_records(set, filters)

    for annotation_record in annotation_records:
        mapping: DataFilesMapping = DataFilesMapping()
        
        records_key = first_or_default(dga_map.keys(), lambda x: annotation_record["was_informed_by"][0] in x)
        if records_key is None:
            log.warning(f"not found: {annotation_record['was_informed_by']}")
            continue
        records = dga_map[records_key]

        set = NMDCCollection.data_object_set
        filters = [
            filter_field_value_matches_any("id", [*annotation_record["has_output"], *records["metaproteome"]["has_output"]]),
            filter_field_value_matches_any("data_object_type", ["Annotation Amino Acid FASTA", "Functional Annotation GFF", "LC-DDA-MS/MS Raw Data"])
        ] 
        data_object_records = get_records(set, filters)
        to_download_list.extend([("data", do) for do in data_object_records])

        index = {do["data_object_type"]: do for do in data_object_records}

        raw_do, gff_do, faa_do = itemgetter(
            "LC-DDA-MS/MS Raw Data",
            "Functional Annotation GFF",
            "Annotation Amino Acid FASTA"
        )(index)

        # build mapping object
        mapping.data_generation_id = colon_to_underscore(records["metaproteome"]["id"])
        mapping.dataset_name = Path(raw_do["name"]).stem
        mapping.raw_file_loc = datafiles_path / raw_do["name"]
        mapping.faa_file_loc = datafiles_path / faa_do["name"]
        mapping.gff_file_loc = datafiles_path / gff_do["name"]
        mapping.dataset_id = colon_to_underscore(raw_do["id"])
        mapping.faa_file_id = colon_to_underscore(faa_do["id"])
        mapping.gff_file_id = colon_to_underscore(gff_do["id"])
        data_files_mapping.append(mapping)


    # get settings files
    # msgf+ param file
    set = NMDCCollection.data_object_set
    filters = [filter_field_value_equals("id", msgf_param_id)] 
    data_object_records = get_records(set, filters)
    msgf_do = data_object_records[0]
    to_download_list.append(("settings", msgf_do))
    
    # masic param file
    set = NMDCCollection.data_object_set
    filters = [filter_field_value_equals("id", masic_param_id)] 
    data_object_records = get_records(set, filters)
    masic_do = data_object_records[0]
    to_download_list.append(("settings", masic_do))

    # contam file
    set = NMDCCollection.data_object_set
    filters = [filter_field_value_equals("id", contam_id)] 
    data_object_records = get_records(set, filters)
    contam_do = data_object_records[0]
    to_download_list.append(("settings", contam_do))

    # build input.json
    metap_input: MetapInput = MetapInput()
    metap_input.mapper_list = data_files_mapping
    metap_input.masic_param_file_filepath = settings_path / masic_do["name"]
    metap_input.msgf_param_file_filepath = settings_path / msgf_do["name"]
    metap_input.kaiko_param_file_filepath = settings_path / "kaiko_defaults.yaml"
    metap_input.contaminant_param_file_filepath = settings_path / contam_do["name"]
    metap_input.masic_parameter_file_id = colon_to_underscore(masic_do["id"])
    metap_input.msgf_parameter_file_id = colon_to_underscore(msgf_do["id"])
    metap_input.contaminant_parameter_file_id = colon_to_underscore(contam_do["id"])
    metap_input.qvalue_threshold = q_value_threshold
    metap_input.study = colon_to_underscore(study_id)
    metap_input.execution_resource = execution_resource
    metap_input.data_url = data_url
    metap_input.metagenome_free = False # for now

    downloaded_files: List[DownloadedFile] = []
    if dl:
        download_files(to_download_list, settings_dir=settings_path, datafiles_dir=datafiles_path, check_md5=md5)

    # write input.json
    with open(output_path / f"{colon_to_underscore(study_id)}_input.json", 'w', encoding="utf-8") as file:
        file.write(metap_input.model_dump_json(indent=4))

    # optionally write downloaded files report
    if dl:
        with open(output_path / "downloaded_files.txt", "w", encoding="utf-8", newline="") as file:
            header = [f.name for f in fields(DownloadedFile)]
            csv_out = csv.DictWriter(file, fieldnames=header)
            csv_out.writeheader()
            csv_out.writerows([asdict(dl) for dl in downloaded_files])


if __name__ == "__main__":
    main()
