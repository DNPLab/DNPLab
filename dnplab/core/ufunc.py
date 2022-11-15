from .data import DNPData
import numpy as _np

__all__ = ["generate_data"]


def generate_data(shape):
    '''Generate random data set with given shape

    Args:
        shape (tuple): shape of data
    
    Returns:
        DNPData: data of given shape
    
    Examples:

        Generate a random 5 x 10 x 20 data set:

        >>> data = generate_data((5, 10, 20))
    '''
    size = 1
    dims = []
    coords = []
    for ix, length in enumerate(shape):
        dims.append("x" + str(ix))
        coords.append(_np.array(range(length)))
        size *= length

    values = _np.random.randn(size)

    return DNPData(values, dims, coords)


def ones(shape, dtype=None):
    values = _np.ones(shape, dtype=dtype)
    coords = []
    dims = []
    for ix in range(len(shape)):
        dims.append(str(ix))
        coords.append(_np.arange(shape[ix]))

    return DNPData(values, dims, coords)


def ones_like(a):
    return ones(a.shape, a.dtype)


def zeros(shape, dtype=None):
    values = _np.zeros(shape, dtype=dtype)
    coords = []
    dims = []
    for ix in range(len(shape)):
        dims.append(str(ix))
        coords.append(_np.arange(shape[ix]))

    return DNPData(values, dims, coords)


def zeros_like(a):
    zeros_ = zeros(a.shape, a.dtype)
    zeros_.dims = a.dims
    zeros_.coords = a.coords
    zeros_.attrs = a.attrs
    return zeros(a.shape, a.dtype)


def randn(shape):
    values = _np.random.randn(*shape)
    coords = []
    dims = []
    for ix in range(len(shape)):
        dims.append(str(ix))
        coords.append(_np.arange(shape[ix]))

    return DNPData(values, dims, coords)


def randn_like(a):
    return randn(a.shape)
