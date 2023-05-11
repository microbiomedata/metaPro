import oauthlib
import requests_oauthlib
import functools

from datetime import datetime
from typing import List
from json import dumps


class NmdcIdSource():
    def __init__(self, client_id=None, client_secret=None, buffer_time_mins=5):
        self._client_id = client_id
        self._client_secret = client_secret
        self._nmdc_mint_url = "https://api.microbiomedata.org/pids/mint"
        self._nmdr_token_url = "https://api.microbiomedata.org/token"
        self._buffer_time_mins = buffer_time_mins
        self._token = None
        self._expires_in_mins = None
        self._init_timestamp = None
        self._client = oauthlib.oauth2.BackendApplicationClient(client_id=self.client_id)
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
                self._token = self._oauth.fetch_token(token_url=self._nmdr_token_url,
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
