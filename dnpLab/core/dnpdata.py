import nddata
import numpy as np


class dnpdata(nddata.nddata_core):

    def __init__(self, values = np.r_[[]], dims = [], coords = [], attrs = {}, error = None, **kwargs):

        super().__init__(values, dims, coords, attrs, error, **kwargs)

    def __repr__(self):
        return 'dnpdata(values = {}, coords = {}, dims = {}, attrs = {})'.format(repr(self._values), repr(self._dims), repr(self._coords), repr(self._attrs))



if __name__ == '__main__':
    x = np.r_[0:10]
    y = np.r_[0:20]
    z = np.r_[0:15]
    data = dnpdata(np.array(range(len(x)*len(y)*len(z))).reshape(len(x),len(y),len(z)), ['x','y','z'], [x, y, z])

