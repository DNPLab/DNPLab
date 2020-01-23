# odnpLab

[![N|Solid](http://www.bridge12.com/wp-content/uploads/2016/10/b12logo.png)](http://www.bridge12.com/)

odnpLab is an open-source, python library for importing and processing ODNP data.

# Features

  - Import ODNP Data from Bruker, Varian, and Kea formats
  - Construct n-dimensional arrays easily
  - Process Data including window, zero-filling, Fourier transforms, etc.

# Installation

  - A description of the installation procedure will be available on release

# Developement 

  - Clone the git repository from https://bitbucket.org/tmaly/odnplab/

# Importing Library without installing

While in development, the code will not be installed in the same location as other python packages. To import the package, you must direct Python to look in the folder containing the package on your computer.

```python
import sys
sys.path.append('path/to/odnpLab/package')
# 'path/to/odnpLab/package' is a path to the folder containing the README.md and setup.py files
```

# Example Script for Bruker Data

```python
# Import standard packages
import numpy as np
from matplotlib.pylab import *

# Append odnpLab Package to Python Path
import sys
sys.path.append('../odnpLab')

# Import odnpLab package
import odnpLab

# Define Path to Data
path = '../data/HanLab_odnpData/2017_01_10_AMS37_DNP_hydrophilic/'

# Import data from given expNum
data =  odnpLab.odnpImport.bruker.importBruker(path,expNum = 20)

# Define data dictionary
dataDict = {}
dataDict['raw'] = data

# Define Dictionary of processing parameters
procDict = {}

dataDict = odnpLab.odnpNMR.window(dataDict,procDict)

figure('data raw')
dataDict['raw'].plot('t',linewidth = 2.0)
xlabel('Time (s)')
ylabel('Signal (a.u.)')

figure('data window')
dataDict['proc'].plot('t',linewidth = 2.0)
xlabel('Time (s)')
ylabel('Signal (a.u.)')
dataDict = odnpLab.odnpNMR.fourierTransform(dataDict,procDict)

# How to access data stored in numpy array and normalize data
dataDict['proc'].data = dataDict['proc'].data / np.max(np.abs(dataDict['proc'].data))

# Autophase routine built into odnpData class
dataDict['proc'].autophase()

# Select data over range of -50 to +50 ppm
dataDict['proc'] = dataDict['proc'].range('t',-50,50)

figure('spectrum')
dataDict['proc'].plot('t',markersize = 8.)
show()
```

# Example Script for Kea Data

```python
# Import standard packages
import numpy as np
from matplotlib.pylab import *

# Append odnpLab path to Python Path
import sys
sys.path.append('G:/My Drive/Exchange/Projects/0055 CPF-NIGM-0055 ODNP System/Software/Python/odnpLab')

# import odnpLab
import odnpLab

# Define data Directory
path = 'G:/My Drive/Exchange/Projects/0055 CPF-NIGM-0055 ODNP System/Software/Python/ODNP Class/ref_data/10mM_TEMPOL_water_11-20-2019'

# Experiment List
exp_num_list = [1,24,33,37,39]

# Convert to Array
exp_num_array = np.array(exp_num_list)

# Convert to power in dBm
power_dBm_array = exp_num_array - 2.
cable_loss = 1.5 # cable power loss in dB
power_dBm_array -= cable_loss

# Convert to Watts
power_W_array = 10.**((power_dBm_array-30.)/10.)
power_W_array[0] = 0 # Correct first value

# Pre-allocate odnpData object
data_power = odnpLab.odnpData()
# Import Kea data and concatenate power dimension
for ix in range(len(exp_num_list)):
    temp_data = odnpLab.odnpImport.kea.importKea(path,num = exp_num_list[ix])
    data_power.add_power(temp_data,power_W_array[ix])

allData = {'raw':data_power}

# Fourier Transform data
allData = odnpLab.odnpNMR.fourierTransform(allData,{})

figure('FT data')
allData['proc'].plot('t')
show()
```
