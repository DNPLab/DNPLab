import sys
sys.path.append('..')
import numpy as np

import dnpLab as dnp

filename = 'test.h5'

ws = dnp.dnpImport.h5.loadh5(filename)

dnp.plot(ws['ft'])
dnp.xlim(20,-20)
dnp.figure()
dnp.plot(ws['proc'], 'o')

dnp.figure
ws['ft'].reorder(['power'])
dnp.imshow(ws['ft'].abs)

dnp.show()

