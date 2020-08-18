# dnpLab - Bringing the Power of Python to DNP-NMR Spectroscopy

dnpLab is a collaboration between:
- [Bridge12 Technologies, Inc.](http://www.bridge12.com/)
- [Han Lab UCSB](https://han.chem.ucsb.edu/)
- [Franck Lab Syracuse University](https://jmfrancklab.github.io/)

dnpLab is an open-source, python library for importing and processing DNP-NMR data.

# Features

  - Import DNP-NMR data from Topspin, VnmrJ, and Prospa formats
  - Construct N-dimensional data objects easily
  - Processing data: including window, zero-filling, Fourier transforms, alignment, etc.

# Requirements

  - Python 3.6 or Later
  - odnpLab requires: numpy, matplotlib, scipy, h5py, PyQt5

To install the required packages, perform:
```console
pip install numpy matplotlib scipy h5py PyQt5
```

# Installation

dnpLab can be installed via pip:

```console
pip install dnpLab
```

# Developement 

  - Clone the git repository from https://bitbucket.org/tmaly/odnplab/

# Example Script for Topspin Data

```python
# Import odnpLab package
import dnpLab as dnp

# Define Path to Data
path = 'path/to/data'

# Import data from given expNum
data =  dnp.dnpImport.topspin.import_topspin(path, expNum = 20)

# Create workspace for data
ws = dnp.create_workspace('raw', data)
ws.copy('raw', 'proc') # copy data to processing buffer ('proc' by default)

# Exponential Apodization with 20 Hz linewidth
ws = odnpLab.odnpNMR.window(ws, {'linewidth': 20.})

# Create figure
dnp.dnpResults.plt.figure()
dnp.dnpResults.plot(ws['raw'])
dnp.dnpResults.plot(ws['proc'])

# Add axis labels
dnp.dnpResults.plt.xlabel('Time (s)')
dnp.dnpResults.plt.ylabel('Signal (a.u.)')

# Perform Fourier Transform with default parameters
ws = odnpLab.odnpNMR.fourierTransform(ws, {})

# Create figure for spectrum
dnp.dnpResults.plt.figure()
dnp.dnpResults.plot(ws['proc'])

# Add axis labels
dnp.dnpResults.plt.xlabel('Chemical Shift (ppm)')
dnp.dnpResults.plt.ylabel('Signal (a.u.)')

# set x-limits from -100 to 100 ppm
dnp.dnpResults.plt.xlim(100, -100)
```
