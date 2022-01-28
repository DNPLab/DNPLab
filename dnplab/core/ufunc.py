from .data import DNPData
import numpy as np

__all__ = ['generate_data']

def generate_data(shape):
    size = 1
    dims = []
    coords = []
    for ix, length in enumerate(shape):
        dims.append('x' + str(ix))
        coords.append(np.array(range(length)))
        size *= length

    values = np.random.randn(size)

    return DNPData(values, dims, coords)


def ones(shape, dtype=None):
    values = np.ones(shape, dtype=dtype)
    coords = []
    dims = []
    for ix in range(len(shape)):
        dims.append(str(ix))
        coords.append(np.arange(shape[ix]))

    return DNPData(values, dims, coords)


def ones_like(a):
    return ones(a.shape, a.dtype)


def zeros(shape, dtype=None):
    values = np.zeros(shape, dtype=dtype)
    coords = []
    dims = []
    for ix in range(len(shape)):
        dims.append(str(ix))
        coords.append(np.arange(shape[ix]))

    return DNPData(values, dims, coords)


def zeros_like(a):
    zeros_ = zeros(a.shape, a.dtype)
    zeros_.dims = a.dims
    zeros_.coords = a.coords
    zeros_.attrs = a.attrs
    return zeros(a.shape, a.dtype)


def randn(shape):
    values = np.random.randn(*shape)
    coords = []
    dims = []
    for ix in range(len(shape)):
        dims.append(str(ix))
        coords.append(np.arange(shape[ix]))

    return DNPData(values, dims, coords)


def randn_like(a):
    return randn(a.shape)
