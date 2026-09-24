import urllib.request
import shutil
import os
from pathlib import Path
from typing import Optional


def download_file(url: str, base_path: str) -> Optional[str]:
    to_return = None
    file_name = os.path.basename(url)
    filepath = Path(base_path) / file_name

    print(f"File downloading {file_name} from {url}")

    with urllib.request.urlopen(url) as response:
        if response.status == 200:
            with open(filepath, 'wb') as fp:
                shutil.copyfileobj(response, fp)
            print(f"File downloaded and saved as {filepath}")
            to_return = str(filepath.absolute())
        else:
            print(f"Failed to download file. Status code: {response.status}")
            to_return = None

    return to_return