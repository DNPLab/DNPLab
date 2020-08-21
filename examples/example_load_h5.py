import sys
sys.path.append('..')
import numpy as np

import dnplab as dnp

filename = 'test.h5'

ws = dnp.dnpImport.h5.loadh5(filename)

dnp.dnpResults.plot(ws['ft'])
dnp.dnpResults.xlim(20,-20)
dnp.dnpResults.figure()
dnp.dnpResults.plot(ws['proc'], 'o')

dnp.dnpResults.figure()
ws['ft'].reorder(['power'])
dnp.dnpResults.imshow(ws['ft'].abs)

dnp.dnpResults.show()

