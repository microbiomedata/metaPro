import oauthlib
import requests_oauthlib
import functools
import requests
import json
import random as r
import string as s

from datetime import datetime
from typing import List
from json import dumps


class NmdcIdSource():
    def __init__(self, client_id=None, client_secret=None, username = None, password = None, buffer_time_mins=5):
        self._client_id = client_id
        self._client_secret = client_secret
        self._username = username
        self._password = password
        self._nmdc_mint_url = "https://api.microbiomedata.org/pids/mint"
        self._nmdc_token_url = "https://api.microbiomedata.org/token"
        self._buffer_time_mins = buffer_time_mins
        self._token = None
        self._expires_in_mins = None
        self._init_timestamp = None
        self._client = oauthlib.oauth2.BackendApplicationClient(client_id=self._client_id)
        self._oauth = requests_oauthlib.OAuth2Session(client=self._client)

    def _save_token(self, token):
        self._token = token
        self._init_timestamp = datetime.now()
        self._expires_in_mins = token['expires']['minutes']

    def _get_extras(self):
        return {
            'client_id': self._client_id,
            'client_secret': self._client_secret,    
        }

    def _check_and_refresh(method):
        @functools.wraps(method)
        def check_and_refresh(self, *args, **kwargs):
            if self._init_timestamp == None or ((datetime.now() - self._init_timestamp).total_seconds() / 60 > (self._expires_in_mins - self._buffer_time_mins)):
                self._token = self._oauth.fetch_token(token_url=self._nmdc_token_url,
                                                      client_id=self._client_id,
                                                      client_secret=self._client_secret)
                self._init_timestamp = datetime.now()
                self._expires_in_mins = self._token['expires']['minutes']
            return method(self, *args, **kwargs)
            
        return check_and_refresh

    @_check_and_refresh
    def get_ids(self, type: str, how_many: int) -> List[str]:
        payload = {
            'schema_class': {
                'id': type,
            },
            'how_many': how_many
        }
        response = self._oauth.post(self._nmdc_mint_url, data=dumps(payload))

        if(response.status_code != 200):
            raise Exception(f"Error getting ids from NMDC: {response.status_code} {response.text}")

        return response.json()
    
    @staticmethod
    def get_token_from_password(username: str, password: str) -> str:
        nmdc_token_url = 'https://api.microbiomedata.org/token'
        headers = {
            'accept': 'application/json',
            'Content-Type': 'application/x-www-form-urlencoded'
        }
        data = {
            'grant_type': 'password',
            'username': username,
            'password': password,
            'scope': '',
            'client_id': '',
            'client_secret': ''
        }

        response = requests.post(nmdc_token_url, data=data, headers=headers)
        if response.status_code == 200:
            return response.json()
        else:
            print(f"Error getting token from NMDC: {response.status_code} {response.text}")
            return None

    @staticmethod
    def get_ids_from_token(type: str, how_many: int, token: str):
            url = "https://api.microbiomedata.org/pids/mint"
            headers = {
                'accept': 'application/json',
                'Authorization': f'Bearer {token}',
                'Content-Type': 'application/json'
            }
            payload = {
                'schema_class': {
                    'id': type,
                },
                'how_many': how_many
            }
            response = requests.post(url, data=payload, headers=headers)

            if(response.status_code != 200):
                raise Exception(f"Error getting ids from NMDC: {response.status_code} {response.text}")

            return response.json()


class NmdcFakeIdSource:
    """
        A Class for generating fake IDs for testing
    """
    def __init__(self, client_id=None, client_secret=None, username = None, password = None):        
        self._client_id = client_id
        self._client_secret = client_secret
        self._username = username
        self._password = password
  
    def get_ids(self, type: str, how_many: int) -> List[str]:
        type_code: str = ""
        shoulder = "1"

        if type == "nmdc:DataObject":
            type_code = "dobj"
        elif type == "nmdc:MetaproteomicsAnalysisActivity":
            type_code = "wfmp"
        else:
            type_code = "unkn"

        ids = [
            f"nmdc:{type_code}-{shoulder}-xxx" + ''.join(r.choices(s.ascii_lowercase + s.digits, k=5))
            for _ 
            in range(how_many)
        ]

        return ids