# dnplab example script to process inversion recovery data (Bruker format)
 
# Import python modules
# If dnplab is installed using pip, the first two lines are not required
import sys
import numpy as np

sys.path.append('..')
import dnplab

# Import data
path = '\path\to\topsin\data\'
folder = 304
data = dnplab.dnpImport.topspin.import_topspin(path + os.sep,folder)

# Create workspace
ws = dnplab.create_workspace()
ws.add('raw', data)
ws.copy('raw', 'proc')

# Process FIDs
dnplab.dnpNMR.remove_offset(ws,{})
dnplab.dnpNMR.window(ws,{'linewidth' : 10})
dnplab.dnpNMR.fourier_transform(ws,{'zero_fill_factor' : 2})
dnplab.dnpNMR.align(ws, {})
dnplab.dnpNMR.autophase(ws,{})

ws.copy('proc', 'ft')
dnplab.dnpNMR.integrate(ws, {'integrate_width' : 100, 'integrate_center' : 0})

# Display NMR spectra
dnplab.dnpResults.plot(ws['ft'].real)
dnplab.dnpResults.xlim([-30,50])
dnplab.dnpResults.plt.xlabel('Chemical Shift [ppm]')
dnplab.dnpResults.plt.ylabel('Signal Amplitude [a.u.]')

# Fit inversion-recovery build-up, display T1 and plot results
dnplab.dnpFit.t1Fit(ws)
print('T1 value (sec) = ' + str(ws['fit'].attrs['t1']))

<<<<<<< HEAD
dnplab.dnpResults.figure()
dnplab.dnpResults.plot(ws['proc'].real, 'o')
dnplab.dnpResults.plot(ws['fit'])
dnplab.dnpResults.plt.xlabel('Time t1 [s]')
dnplab.dnpResults.plt.ylabel('Signal Amplitude [a.u.]')
dnplab.dnpResults.show()
=======
dnp.dnpResults.figure()
dnp.dnpResults.plot(ws['proc'].real, 'o')
dnp.dnpResults.plot(ws['fit'])
dnp.dnpResults.plt.xlabel('Time t1 [s]')
dnp.dnpResults.plt.ylabel('Signal Amplitude [a.u.]')
dnp.dnpResults.show()
>>>>>>> develop
