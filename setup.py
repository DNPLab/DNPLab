import setuptools
from distutils.core import setup

with open("README.md", "r") as f:
    long_description = f.read()

with open("dnplab/version.py", "r") as f:
    # Define __version__
    exec(f.read())

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
        "Examples": "http://dnplab.net/auto_examples/index.html",
        "Source Code": "https://github.com/DNPLab/dnplab",
        "Web App": "http://odnplab.com/",
    },
    keywords=["ODNP", "DNP", "NMR"],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.24",
        "scipy>=1.10",
        "matplotlib>=3.7",
        "h5py>=3.8",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    package_data={"dnplab", ["config/dnplab.cfg"]},
)
