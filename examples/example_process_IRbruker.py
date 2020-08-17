# dnpLab example script to process inversion recovery data (Bruker format)
 
# Import python modules
# If dnpLab is installed using pip, the first two lines are not required
import sys
sys.path.append('..')

import numpy as np
import dnpLab as dnp

# Import data
path = '\path\to\topsin\data\'
folder = 304
data = dnp.dnpImport.topspin.import_topspin(path + os.sep,folder)

# Create workspace
ws = dnp.create_workspace()
ws.add('raw', data)
ws.copy('raw', 'proc')

# Process FIDs
dnp.dnpNMR.remove_offset(ws,{})
dnp.dnpNMR.window(ws,{'linewidth' : 10})
dnp.dnpNMR.fourier_transform(ws,{'zero_fill_factor' : 2})
dnp.dnpNMR.align(ws, {})
dnp.dnpNMR.autophase(ws,{})

ws.copy('proc', 'ft')
dnp.dnpNMR.integrate(ws, {'integrate_width' : 100, 'integrate_center' : 0})

# Display NMR spectra
dnp.dnpResults.plot(ws['ft'].real)
dnp.dnpResults.xlim([-30,50])
dnp.dnpResults.plt.xlabel('Chemical Shift [ppm]')
dnp.dnpResults.plt.ylabel('Signal Amplitude [a.u.]')

# Fit inversion-recovery build-up, display T1 and plot results
dnp.dnpFit.t1Fit(ws)
print('T1 value (sec) = ' + str(ws['fit'].attrs['t1']))

dnp.dnpResults.figure()
dnp.dnpResults.plot(ws['proc'].real, 'o')
dnp.dnpResults.plot(ws['fit'])
dnp.dnpResults.plt.xlabel('Time t1 [s]')
dnp.dnpResults.plt.ylabel('Signal Amplitude [a.u.]')
dnp.dnpResults.show()