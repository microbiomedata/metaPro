import fastjsonschema
import json
from jsonschema import validate
import jsonschema
import requests



def check_data_object(d):
    """Validate the metadata urls.

    Args:
        d ([type]): [description]

    Returns:
        [type]: [description]
    """

    rv = requests.head(
        d["url"],
        allow_redirects=True,
        verify=False,
        timeout=5,
        headers={"Accept-Encoding": "gzip;q=0"},
    )

    if not rv.status_code == 200:
        return {
            "error": {"status_code": rv.status_code, "details": "not OK"},
            "id": d["id"],
        }

    if d["file_size_bytes"] != int(rv.headers["Content-Length"]):
        return {
            "error": {
                "details": "file size different than reported",
                "file_size_actual": rv.headers["Content-Length"],
                "file_size_reported": d["file_size_bytes"],
            },
            "id": d["id"],
        }

    return {"result": "OK", "id": d["id"]}


def using_fastjsonschema(schema, json_file):
    with open(schema) as json_schema:
        # load schema.
        nmdc_schema = json.load(json_schema)
        try:
            schema_validator = fastjsonschema.compile(nmdc_schema)
        except fastjsonschema.JsonSchemaDefinitionException as bad_def:
            print({"bad_def": bad_def})
        # validate data.
        with open(json_file) as json_data:
            try:
                data_objects = json.load(json_data)
            except ValueError:
                print("Invalid json data.")
                raise

            try:
                schema_validator(data_objects)
            except fastjsonschema.JsonSchemaException as bad_data:
                print({"bad_data": bad_data})
                raise

def using_jsonschema(schema, json_file):
    with open(schema) as json_schema, open(json_file) as json_data:
        nmdc_schema = json.load(json_schema)
        try:
            data_objects = json.load(json_data)
        except ValueError:
            print("Can't load json file.")
            raise
        try:
            # for obj in data_objects:

            validate(instance=data_objects, schema=nmdc_schema)
            # break
        except jsonschema.exceptions.ValidationError as err:
            print(err)
            err = "JSON data is InValid"
            return False, err

def validate_file(filename):
    capture_error = []
    counter = 0
    with open(filename) as f:
        data = json.load(f)
        for key, val in data.items():
            for obj in val:
                counter += 1
                print(check_data_object(obj))
            print(f"{counter} number of objects parsed.")

if __name__ == "__main__":
    # activity_json_file = "storage/results/stegen/stegen_MetaProteomicAnalysis_activity.json"
    # data_obj_json_file = "storage/results/stegen/stegen_emsl_analysis_data_objects.json"

    activity_json_file = "mpdocs/temp/stegen_MetaProteomicAnalysis_activity.json"
    data_obj_json_file = "mpdocs/temp/stegen_emsl_analysis_data_objects.json"

    # schema = "./metaPro/docs/nmdc.schema.json"
    schema = "docs/nmdc.schema.json"
    evaluate_file = data_obj_json_file

    # nmdc schema validation
    print("activity_json schema validation")
    using_fastjsonschema(schema, activity_json_file)
    print("data_obj_json schema validation")
    using_fastjsonschema(schema, data_obj_json_file)
     # # using_jsonschema(schema, evaluate_file)
    print('<<>>' * 50)
    # urls correctness validation
    print("data_obj_json Url correctness")
    validate_file(data_obj_json_file)

