"""
Download dataset from LOGS

"""

import os
import zipfile
from LOGS import LOGS
from LOGS.Entities import DatasetRequestParameter
from . import load


def download(
    names: str | list[str], url: str, apiKey: str | None = None, verify: bool = False
):
    """
    Download dataset from LOGS

    """
    if isinstance(names, str):
        names = [names]

    request = DatasetRequestParameter(names=names)

    logs = LOGS(url=url, apiKey=apiKey, verify=verify)

    datasets = logs.datasets(request)
    if not datasets:
        raise ValueError("No datasets found with the given parameters.")

    if datasets.count != len(names):
        raise ValueError(
            "The number of datasets found does not match the number of names provided."
        )

    for d in datasets:
        if not d.claimed:
            print("Warning: Please claim the dataset")
            break

    path = os.path.join(os.getcwd(), "data")
    zipfilename = "LOGS.zip"

    if "data" not in os.listdir(os.getcwd()):
        os.mkdir(path)

    datasets.download(path, zipfilename, overwrite=True)

    zipFile = os.path.join(path, zipfilename)

    with zipfile.ZipFile(zipFile, "r") as zip_ref:
        zip_ref.extractall(path)

    data_format = None
    for file in os.listdir(path):
        path_exten = os.path.splitext(file)[1]
        if path_exten == ".zip":  # ignore zip files
            continue
        try:
            data_format = load.autodetect(file)
            break
        except TypeError:
            continue

    file_list = [
        os.path.join(path, file) for file in os.listdir(path) if path_exten in file
    ]

    return file_list, data_format
