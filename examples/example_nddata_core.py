import sys
sys.path.append('..')

<<<<<<< HEAD
from dnpLab.core.nddata import nddata_core
import dnplab
=======
from dnplab.core.nddata import nddata_core
import dnplab as dnp
>>>>>>> develop

import numpy as np

x = np.r_[1:100]
values = np.random.randn(len(x))

data = nddata_core(values, ['x'], [x])
print(data)

dnplab.dnpResults.plot(data)
dnplab.dnpResults.show()


