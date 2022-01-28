[![http://dnplab.net](https://img.shields.io/pypi/v/dnplab)](https://pypi.org/project/dnplab/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/dnplab)](https://www.python.org/downloads/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/dnplab.svg?label=Pypi%20downloads)](https://pypi.org/project/dnplab/)
[![Downloads](https://pepy.tech/badge/dnplab/month)](https://pepy.tech/project/dnplab)

[![DNPLab Logo](docs/source/_static/images/dnpLabLogo.png)](http://dnplab.net)
# DNPLab - Bringing the Power of Python to DNP-NMR Spectroscopy

DNPLab is a collaboration between:
- [Bridge12 Technologies, Inc.](http://www.bridge12.com/)
- [Han Lab UCSB](https://han.chem.ucsb.edu/)
- [Franck Lab Syracuse University](https://jmfrancklab.github.io/)

Authors:
Timothy Keller, Thomas Casey, Yanxian Lin, John Franck, Thorsten Maly, Songi Han

DNPLab is an open-source, python library for importing and processing DNP-NMR data.

# Features

  - Import DNP-NMR data from Topspin, VnmrJ, and Prospa formats
  - Construct N-dimensional data objects easily
  - Process data using apodization, zero-filling, Fourier transformations, alignment, etc.
  - Analyze data using the hydration module

# Requirements

  - Python 3.6 or Later
  - DNPLab requires: numpy, matplotlib, scipy, h5py, PyQt5

To install the required packages, perform:
```console
pip install numpy matplotlib scipy h5py PyQt5
```

# Installation

DNPLab can be installed via pip:

```console
pip install dnplab
```

# Documenation

Find the online documentation at: http://docs.dnplab.net

# Developement 

  - Clone the git repository from https://github.com/DNPLab/dnpLab