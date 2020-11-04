from distutils.core import setup

import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

with open("dnplab/version.py", "r") as f:
    # Define __version__
    __version__ = eval(f.read().split("=")[-1].strip())

setup(
    name="dnplab",
    packages=setuptools.find_packages(),
    version=__version__,
    license="MIT",
    description="DNPLab - Bringing the Power of Python to DNP-NMR Spectroscopy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="DNPLab Team",
    url="http://dnpLab.net",
    download_url="",
    project_urls={
        "Documentation": "http://docs.dnpLab.net",
        "Examples": "http://examples.dnplab.net/",
        "Source Code": "https://github.com/DNPLab/dnplab",
        "Web App": "http://odnplab.com/",
    },
    keywords=["ODNP", "DNP", "NMR"],
    python_requires=">=3.6",
    install_requires=[
        "numpy>=1.19",
        "scipy>=1.5",
        "matplotlib>=3.3",
        "h5py>=2.10",
        "PyQt5==5.13",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    entry_points=dict(
        console_scripts=[
            "hydrationGUI=dnplab.hydrationGUI:main_func",
        ]
    ),
)
