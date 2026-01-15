import abc
import argparse
import csv
import json
import sys
import sqlite3
import time
from dataclasses import dataclass, fields, asdict, astuple
from itertools import chain
from pathlib import Path
import urllib.request
from typing import List, Union


class NMDCCollection:
    biosample_set: str = "biosample_set"
    calibration_set: str = "calibration_set"
    collecting_biosamples_from_site_set: str = "collecting_biosamples_from_site_set"
    configuration_set: str = "configuration_set"
    data_generation_set: str = "data_generation_set"
    data_object_set: str = "data_object_set"
    field_research_site_set: str = "field_research_site_set"
    functional_annotation_agg: str = "functional_annotation_agg"
    functional_annotation_set: str = "functional_annotation_set"
    genome_feature_set: str = "genome_feature_set"
    instrument_set: str = "instrument_set"
    manifest_set: str = "manifest_set"
    material_processing_set: str = "material_processing_set"
    processed_sample_set: str = "processed_sample_set"
    storage_process_set: str = "storage_process_set"
    study_set: str = "study_set"
    workflow_execution_set: str = "workflow_execution_set"


class FilterExpression(abc.ABC):
    def __init__(self, field: str, value: str):
        self.field = field
        self.value = value

    def get_expression(self) -> str:
        pass


class FilterExpressionEqual(FilterExpression):
    def __init__(self, field: str, value: str):
        super().__init__(field, value)

    def get_expression(self) -> str:
        return f'"{self.field}": "{self.value}"'


class FilterExpressionIn(FilterExpression):
    def __init__(self, field: str, values: List[str], spec: str):
        self.spec = spec
        super().__init__(field, values)

    def get_expression(self) -> str:
        values_str = ', '.join([f'"{value}"' for value in self.value])
        return f'"{self.field}": {{ "{self.spec}": [{values_str}] }}'


class Filter:
    def __init__(self, filter_expressions: List[FilterExpression] = []):
        self.filter_expressions = filter_expressions

    def get_filter_query(self):
        expression_str = ', '.join([f'{value.get_expression()}' for value in self.filter_expressions])
        return '{' + expression_str + '}'
    
    def with_expression(self, expression: FilterExpression):
        self.filter_expressions.append(expression)
        return self
    
    def with_field_value_equals(self, field: str, to_compare: str):
        self.filter_expressions.append(FilterExpressionEqual(field, to_compare))
        return self
    
    def with_field_value_matches_all(self, field: str, to_compare_values: List[str]):
        self.filter_expressions.append(FilterExpressionIn(field, to_compare_values, "$all"))
        return self
    
    def with_field_value_matches_any(self, field: str, to_compare_values: List[str]):
        self.filter_expressions.append(FilterExpressionIn(field, to_compare_values, "$in"))
        return self
    
    @staticmethod
    def field_value_equals(field: str, to_compare: str):
        return FilterExpressionEqual(field, to_compare)

    @staticmethod
    def field_value_matches_all(field: str, to_compare_values: List[str]):
        return FilterExpressionIn(field, to_compare_values, "$all")

    @staticmethod
    def field_value_matches_any(field: str, to_compare_values: List[str]):
        return FilterExpressionIn(field, to_compare_values, "$in")


def _get(url, filter_query: str, max_page_size: int = 20, next_page_token: str = None):
    base_url = url
    params = {}

    params["filter"] = filter_query
    params["max_page_size"] = max_page_size
    if next_page_token:
        params["page_token"] = next_page_token

    headers = {
        'accept': 'application/json'
        }
    
    full_url = f"{base_url}?{urllib.parse.urlencode(params)}"
    
    request = urllib.request.Request(full_url, headers=headers, method='GET')

    try:
        with urllib.request.urlopen(request) as response:
            response_body = response.read().decode('utf-8')
            if response.status == 200:
                return json.loads(response_body)
            else:
                print(f"Failed to fetch data for '{filter_query}' -- response code: {response.status} -- response: {response.response_body}")
                return None
    except urllib.error.HTTPError as e:
        print(f"Failed to fetch data for '{filter.get_filter_query()}' -- response code: {e.code} -- response: {e.reason}")
    
    return None


def nmdc_get(url, filter: Filter, max_page_size: int = 60):
    results = []
    next_page_token = None

    while True:
        response = _get(url, filter.get_filter_query(), max_page_size, next_page_token)
        if response:
            result = response.get("resources", [])
            
            if len(result) == 0:
                break

            results.extend(result)
            next_page_token = response.get("next_page_token")
            if not next_page_token:
                break
        else:
            break

    return results


## helpers
def filter_field_value_equals(field: str, to_compare: str):
    return Filter.field_value_equals(field, to_compare)


def filter_field_value_matches_all(field: str, to_compare_values: List[str]):
    return Filter.field_value_matches_all(field, to_compare_values)


def filter_field_value_matches_any(field: str, to_compare_values: List[str]):
    return Filter.field_value_matches_any(field, to_compare_values)


## main access function
def get_records(collection: str, filter_on: Union[Filter|List[FilterExpression]]):
    base_url = f"https://api-backup.microbiomedata.org/nmdcschema/{collection}"
    filter = None
    
    if isinstance(filter_on, Filter):
        filter = filter_on    
    elif isinstance(filter_on, list) and all(isinstance(i, FilterExpression) for i in filter_on):
        filter = Filter(filter_on)
    else:
        raise ValueError("filter_on must be a Filter or a list of FilterExpression")
    
    return nmdc_get(base_url, filter)


parser = argparse.ArgumentParser(
    description="Get metaproteomics ID table from the NMDC database",
    prog='Get Metaproteomics Output Summary',
    epilog='Outputs a summary table of metaproteomics data objects in the NMDC database. This script has no external dependencies.'
)

g = parser.add_mutually_exclusive_group(required=False)

g.add_argument("-t", "--tsv",
               dest="format", action="store_const", const="tsv",
               help="Output as TSV file (default)")
g.add_argument("-s", "--sqlite",
               dest="format", action="store_const", const="sql",
               help="Output as SQLite database")

parser.set_defaults(format="tsv")  # default when neither flag is supplied

@dataclass
class Row:
    study_id: str = ""
    biosample_or_processed_sample_id: str = ""
    data_generation_id: str = ""
    workflow_execution_id: str = ""
    data_object_id: str = ""


def main():
    try:
        args = parser.parse_args()
    except:
        print("Error parsing arguments", file=sys.stderr)
        sys.exit(1)

    format = args.format

    # get all workflow execution records (WFE) for metaproteomics
    filter_mpa = filter_field_value_equals("type", "nmdc:MetaproteomicsAnalysis")
    wfe_records = get_records("workflow_execution_set", [ filter_mpa ])
    wfe_records_map = {record["id"]: record for record in wfe_records}
    
    # Convenience map for data object id -> workflow execution id
    do_id_to_wfe_id_map = {}
    for wfe_record in wfe_records:
        for do_id in wfe_record["has_output"]:
            do_id_to_wfe_id_map[do_id] = wfe_record["id"]

    # get all data generation records associated with WFE records
    filter_dg = filter_field_value_matches_any("id", [record["was_informed_by"][0] for record in wfe_records])
    dg_records = get_records("data_generation_set", [ filter_dg ])
    dg_records_map = {record["id"]: record for record in dg_records}

    # Get all data objects from all WFE outputs
    output_dos = chain.from_iterable([mpa_records["has_output"] for mpa_records in wfe_records])
    filter_do = filter_field_value_matches_any("id", list(output_dos))
    do_records = get_records("data_object_set", [ filter_do ])

    # build rows by working back from data object
    row_list = []
    for do_record in do_records:
        row = Row()
        row.data_object_id = do_record["id"]

        # DO may not always have was_generated_by field
        if "was_generated_by" in do_record:
            row.workflow_execution_id = do_record["was_generated_by"]
        else:
            row.workflow_execution_id = do_id_to_wfe_id_map[do_record["id"]]
        
        row.data_generation_id = wfe_records_map[row.workflow_execution_id]["was_informed_by"][0]
        row.biosample_or_processed_sample_id = dg_records_map[row.data_generation_id]["has_input"][0]    
        row.study_id = dg_records_map[row.data_generation_id]["associated_studies"][0]
        
        row_list.append(row)

    # write output 
    if format == "tsv":
        with open(f"metaproteomics_output_summary.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=[field.name for field in fields(Row)])
            writer.writeheader()
            for row in row_list:
                writer.writerow(asdict(row))

    elif format == "sql":
        db_path = Path(f"metaproteomics_output_summary_{int(time.time() * 1000)}.db")
        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row
        conn.execute('''
            CREATE TABLE IF NOT EXISTS metaproteomics_ids (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                study_id TEXT NOT NULL,
                biosample_or_processed_sample_id TEXT NOT NULL,
                data_generation_id TEXT NOT NULL,
                workflow_execution_id TEXT NOT NULL,
                data_object_id TEXT NOT NULL
            )
        ''')

        rows = [astuple(row) for row in row_list]
        conn.executemany('''
            INSERT INTO metaproteomics_ids (
                study_id,
                biosample_or_processed_sample_id,
                data_generation_id,
                workflow_execution_id,
                data_object_id) Values (?, ?, ?, ?, ?)
        ''', rows)
        conn.commit()
        conn.close()


if __name__ == "__main__":
    main()