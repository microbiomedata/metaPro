import requests
import os
from pathlib import Path
from typing import Optional


def download_file(url, base_path: Path) -> Optional[str]:
    file_name = os.path.basename(url)
    filepath = base_path / file_name

    response = requests.get(url)

    if response.status_code == 200:
        with open(filepath, 'wb') as fp:
            fp.write(response.content)
        print(f"File downloaded and saved as {filepath}")
        return str(filepath.absolute())
    else:
        print(f"Failed to download file. Status code: {response.status_code}")
        return None