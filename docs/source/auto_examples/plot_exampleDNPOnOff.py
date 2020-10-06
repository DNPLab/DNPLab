"""
"This" is my example-script
===========================

This example doesn't do much, it just makes a simple plot
"""

import numpy as np
import dnplab as dnp


########## OFF Signal (P = 0 W) ##########

data = dnp.dnpImport.prospa.import_prospa('../data/prospa/10mM_TEMPO_Water/1Pulse_20200929/35/data.1d')

# data = dnp.dnpImport.prospa.import_prospa('data/1Pulse_20200929/35/data.1d')

ws = dnp.create_workspace()
ws.add('rawOff',data)
ws.copy('rawOff','proc')

dnp.dnpNMR.remove_offset(ws,{})
dnp.dnpNMR.window(ws,{'linewidth' : 15})
dnp.dnpNMR.fourier_transform(ws,{'zero_fill_factor' : 2})

ws.copy('proc','offSignal')


########## ON Signal (P = 0 W) ##########
data = dnp.dnpImport.prospa.import_prospa('../data/prospa/10mM_TEMPO_Water/1Pulse_20200929/51/data.1d')

# data = dnp.dnpImport.prospa.import_prospa('data/1Pulse_20200929/51/data.1d')

ws.add('rawOn',data)
ws.copy('rawOn')

dnp.dnpNMR.remove_offset(ws,{})
dnp.dnpNMR.window(ws,{'linewidth' : 15})
dnp.dnpNMR.fourier_transform(ws,{'zero_fill_factor' : 2})

ws.copy('proc','onSignal')


# ########## Plot Spectra ##########
dnp.dnpResults.figure()
dnp.dnpResults.plot(ws['offSignal'].real * 10 - 100, label = 'Off Signal x 10')
dnp.dnpResults.plot(ws['onSignal'].real, label = 'On Signal')
dnp.dnpResults.xlim([30,-30])
dnp.dnpResults.plt.xlabel('Chemical Shift [ppm]')
dnp.dnpResults.plt.ylabel('Signal Amplitude [a.u.]')
dnp.dnpResults.plt.title('DNP On/Off Signal, 10 mM TEMPO in Water')
dnp.dnpResults.plt.legend()
dnp.dnpResults.plt.grid(True)
dnp.dnpResults.show()