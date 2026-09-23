import abc
import json
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


def filter_field_value_equals(field: str, to_compare: str):
    return Filter.field_value_equals(field, to_compare)


def filter_field_value_matches_all(field: str, to_compare_values: List[str]):
    return Filter.field_value_matches_all(field, to_compare_values)


def filter_field_value_matches_any(field: str, to_compare_values: List[str]):
    return Filter.field_value_matches_any(field, to_compare_values)


def get_records(collection: str, filter_on: Union[Filter|List[FilterExpression]]):
    base_url = f"https://api.microbiomedata.org/nmdcschema/{collection}"
    filter = None
    
    if isinstance(filter_on, Filter):
        filter = filter_on    
    elif isinstance(filter_on, list) and all(isinstance(i, FilterExpression) for i in filter_on):
        filter = Filter(filter_on)
    else:
        raise ValueError("filter_on must be a Filter or a list of FilterExpression")
    
    return nmdc_get(base_url, filter)