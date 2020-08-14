import sys
sys.path.append('..')
import numpy as np

import dnpLab as dnp

path = '..topspin'
folder = 304

data = dnp.dnpImport.topspin.import_topspin(path + os.sep,folder)
ws = dnp.create_workspace()
ws.add('raw', data)
ws.copy('raw', 'proc')

dnp.dnpNMR.remove_offset(ws,{})
dnp.dnpNMR.window(ws,{'linewidth' : 10})
dnp.dnpNMR.fourier_transform(ws,{'zero_fill_factor' : 2})
dnp.dnpNMR.align(ws, {})
dnp.dnpNMR.autophase(ws,{})

ws.copy('proc', 'ft')
dnp.dnpNMR.integrate(ws, {'integrate_width' : 100, 'integrate_center' : 0})

dnp.dnpFit.t1Fit(ws)
print('T1 value (sec) = ' + str(ws['fit'].attrs['t1']))

dnp.dnpResults.plot(ws['ft'])
dnp.dnpResults.xlim([-30,50])
dnp.dnpResults.figure()
dnp.dnpResults.plot(ws['proc'], 'o')
dnp.dnpResults.plot(ws['fit'])
dnp.dnpResults.show()
