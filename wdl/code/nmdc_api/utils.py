import requests
import os
from pathlib import Path
from typing import Optional


def download_file(url, base_path: str) -> Optional[str]:
    to_return = None
    file_name = os.path.basename(url)
    filepath = Path(base_path) / file_name

    print(f"File downloading {file_name} from {url}")

    response = requests.get(url)

    if response.status_code == 200:
        with open(filepath, 'wb') as fp:
            fp.write(response.content)
        print(f"File downloaded and saved as {filepath}")
        to_return = str(filepath.absolute())
    else:
        print(f"Failed to download file. Status code: {response.status_code}")
        to_return = None

    return to_return