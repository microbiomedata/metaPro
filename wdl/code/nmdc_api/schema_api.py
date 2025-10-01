import requests
import os
import abc
from sys import platform
from pathlib import Path
from typing import List, Dict, Optional
from enum import Enum


DEBUG = False


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


def _get(url, filter: Filter, max_page_size: int = 20, next_page_token: str = None):
    base_url = url
    params = {}

    params["filter"] = filter.get_filter_query()
    params["max_page_size"] = max_page_size
    if next_page_token:
        params["page_token"] = next_page_token

    headers = {
        'accept': 'application/json'
        }
    
    response = requests.get(base_url, params=params, headers=headers)
    
    if DEBUG:
        print(response.request.url)
    
    if response.status_code == 200:
        return response.json()  
    else:
        print(f"Failed to fetch data for '{filter.get_filter_query()}' -- response code: {response.status_code} -- response: {response.text}")
        return None

def nmdc_get(url, filter: Filter, max_page_size: int = 60):
    results = []
    next_page_token = None

    while True:
        response = _get(url, filter, max_page_size, next_page_token)
        if response:
            result = response.get("resources", [])
            if DEBUG:
                print(f"rec'd {len(result)} results")
            
            if len(result) == 0:
                break

            results.extend(result)
            next_page_token = response.get("next_page_token")
            if DEBUG:
                print(f"next page token: {next_page_token}")
            if not next_page_token:
                if DEBUG:
                    print("no more pages")
                break
        else:
            break

    return results


### API

def filter_field_value_equals(field: str, to_compare: str):
    return Filter.field_value_equals(field, to_compare)

def filter_field_value_matches_all(field: str, to_compare_values: List[str]):
    return Filter.field_value_matches_all(field, to_compare_values)

def filter_field_value_matches_any(field: str, to_compare_values: List[str]):
    return Filter.field_value_matches_any(field, to_compare_values)

def get_records(collection: str, filter_on: List[FilterExpression]):
    base_url = f"https://api.microbiomedata.org/nmdcschema/{collection}"
    if DEBUG:
        print(f"getting records from {base_url}")
                 
    filter = Filter(filter_on)
    return nmdc_get(base_url, filter)
