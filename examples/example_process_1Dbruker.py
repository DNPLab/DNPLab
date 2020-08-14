import sys
sys.path.append('..')
import numpy as np

import dnpLab as dnp

path = '..topspin'
folder = 20

data = dnp.dnpImport.topspin.import_topspin(path  + os.sep,folder)
ws = dnp.create_workspace()
ws.add('raw', data)
ws.copy('raw', 'proc')

dnp.dnpNMR.remove_offset(ws,{})
dnp.dnpNMR.window(ws,{'linewidth' : 10})
dnp.dnpNMR.fourier_transform(ws,{'zero_fill_factor' : 2})
dnp.dnpNMR.autophase(ws,{})

dnp.dnpResults.figure()
dnp.dnpResults.plot(ws['proc'])
dnp.dnpResults.xlim([-35,50])
dnp.dnpResults.show()
