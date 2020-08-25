import sys
sys.path.append('..')

from dnplab.core.nddata import nddata_core
import dnplab as dnp

import numpy as np

x = np.r_[1:100]
values = np.random.randn(len(x))

data = nddata_core(values, ['x'], [x])
print(data)

dnp.dnpResults.plot(data)
dnp.dnpResults.show()


