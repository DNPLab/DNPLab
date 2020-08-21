# dnpLab example script to import 1D Bruker data, process FID and display NMR spectrum

# Import python modules
# If dnplab is installed using pip, the first two lines are not required
import sys
sys.path.append('..')

import numpy as np
import dnplab as dnp

# Import data
path = '\path\to\topsin\data\'
folder = 20
data = dnp.dnpImport.topspin.import_topspin(path,folder)

# Create workspace
ws = dnp.create_workspace()
ws.add('raw', data)
ws.copy('raw', 'proc')

# Process FID and display spectrum
dnp.dnpNMR.remove_offset(ws,{})
dnp.dnpNMR.window(ws,{'linewidth' : 10})
dnp.dnpNMR.fourier_transform(ws,{'zero_fill_factor' : 2})
dnp.dnpNMR.autophase(ws,{})

dnp.dnpResults.figure()
dnp.dnpResults.plot(ws['proc'].real)
dnp.dnpResults.xlim([-35,50])
dnp.dnpResults.plt.xlabel('Chemical Shift [ppm]')
dnp.dnpResults.plt.ylabel('Signal Amplitude [a.u.]')
dnp.dnpResults.show()
