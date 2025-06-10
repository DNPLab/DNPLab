# %%
"""
.. _plot_03_download_datasets_from_logs:

===============================================
Downloading datasets from LOGS
===============================================

This example demonstrates how to use DNPLab to download datasets from LOGS.


"""
# %%
# ===================================================
# How to download datasets from LOGS
# ===================================================
# This function requires the 'LOGS-py' API to be installed.
#
# To get started, first, setup the python environment:

import dnplab as dnp

# %%
# Lets declare the files you want to download.
# The format can be a string or a list of strings.

files = ['Test Dataset 1', 'Test Dataset 2']

# %%
# U=You might need to provide the URL and API key for LOGS.
url = "https://logs.dnplab.org"
apiKey = 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

# %%
# Next, can download the dataset using the `download` function from the `logs` module.
file_list, data_format = dnp.logs.download(files, url=url, apiKey=apiKey)

# The files will be downloaded to the current working directory in a folder called 'data'.

# Now you can load the files using the `load` function
data = dnp.load(file_list, dim = 't2', coords = [0, 1])
