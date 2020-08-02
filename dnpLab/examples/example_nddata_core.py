import sys
sys.path.append('../..')

from dnpLab.core.nddata import nddata_core
import dnpLab as dnp

import numpy as np

x = np.r_[1:100]
values = np.random.randn(len(x))

data = nddata_core(values, ['x'], [x])

dnp.plot(data)
dnp.show()
print(data)


