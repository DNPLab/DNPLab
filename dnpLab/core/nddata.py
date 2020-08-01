from __future__ import division
import operator
import numpy as np
import warnings
from copy import deepcopy
from collections import OrderedDict
#import nddata_coord
from . import nddata_coord

_numerical_types = (np.ndarray, int, float, complex)

_nddata_core_version = '1.0'

class nddata_core(object):
    '''nddata class
    '''

    __array_priority__ = 1000 # radd, rsub, ... should return nddata object

    def __init__(self, values = np.r_[[]], dims = [], coords = [], attrs = {}, error = None, **kwargs):

        self._nddata_core_version = _nddata_core_version

        # verify values are numpy array
        if isinstance(values, np.ndarray):
            self._values = values
        else:
            raise TypeError('values must be type "numpy.ndarray" not %s'%str(type(values)))
        
        self._coords = nddata_coord.nddata_coord_collection(dims, coords)

        if isinstance(attrs, dict):
            self._attrs = attrs
        else:
            raise TypeError('attrs must be type "dict" not %s'%str(type(attrs)))

        if isinstance(error, np.ndarray) or (error == None):
            self._error = error
        else:
            raise TypeError('error must be type "numpy.ndarray" or "None" not %s'%str(type(error)))

        if 'proc_attrs' in kwargs:
            proc_attrs = kwargs['proc_attrs']
            if isinstance(proc_attrs, list):
                self._proc_attrs = proc_attrs
            else:
                raise TypeError

        if not self._self_consistent():
            warnings.warn('Dimensions not consistent')

    @property
    def __version__(self):
        return self._nddata_core_version

    def _check_dims(self, dims):
        '''Check that list is a list of strings with no duplicates
        
        Args:
            dims list to verify as list of strings

        Returns:
            bool: True if arguments is list of strings with no duplicates, False otherwise
        '''

        # test that all members are strings
        all_strings = all([isinstance(dims[x], str) for x in range(len(dims))])

        # test if any duplicates exist
        any_duplicates = len(dims) == len(set(dims))

        return all_strings and any_duplicates

    def _check_coords(self, coords):
        '''Check that coords is list of 1d numpy arrays
        '''

        for coord in coords:
            if isinstance(coord, np.ndarray):
                if len(coord.shape) == 1:
                    pass
                else:
                    return False
            else:
                return False

        return True

    def _check_error(self, error):
        '''
        '''

        check_type = isinstance(error, np.ndarray)

        if check_type:
            check_size = error.size == self._values.size
        else:
            check_size = False

        return check_type and check_size

    def _self_consistent(self):
        '''Check if dimensions for values, dims, and coords are consistent

        Returns:
            bool: True if consistent, False otherwise
        '''

        if self._values.size == 0:
            coords_check = len(self._coords) == 0
        else:
            coords_check = list(self._values.shape) == list(self.coords.shape)

        dims_check = len(self.values.shape) == len(self.dims)

        return coords_check and dims_check

    def _attrs_valid():
        '''Verify attrs are valid
        '''

        for key in self._attrs:
            if not isinstance(self._attrs, (list, np.ndarray, int, float, str)):
                return False

        return True

    def __getitem__(self, args):
        '''
        '''
        a = self.copy()
        if len(args) % 2 == 1:
            raise ValueError('Cannot index with odd number of arguments')

        index_dims = args[0::2]
        for dim in index_dims:
            if dim not in a.dims:
                raise ValueError('dim not in dims')

        index_slice = args[1::2]

        # check slices
        for slice_ in index_slice:
            # type must be slice or tuple
            if not isinstance(slice_, (slice, tuple, float, int)):
                raise ValueError('Invalid slice type')
            # if tuple, length must be two: (start, stop)
            if isinstance(slice_, tuple) and not len(slice_) in (1, 2):
                raise ValueError('tuple index must have one or two values')

        # convert tuple to slice
        updated_index_slice = []
        for dim, slice_ in zip(index_dims, index_slice):
            if isinstance(slice_, tuple):
                index = a.index(dim)
                if len(slice_) == 1:
                    start = np.argmin(np.abs(slice_[0] - a.get_coord(dim)))
                    updated_index_slice.append(slice(start,start+1))
                else:
                    start = np.argmin(np.abs(slice_[0] - a.get_coord(dim)))
                    stop = np.argmin(np.abs(slice_[1] - a.get_coord(dim)))
                    if start == stop:
                        stop = start + 1
                    if stop < start:
                        start, stop = stop, start
                    updated_index_slice.append(slice(start, stop))
            elif isinstance(slice_, int):
                start = slice_
                if slice_ != -1:
                    updated_index_slice.append(slice(start, start+1))
                else:
                    updated_index_slice.append(slice(slice_, None))
            elif isinstance(slice_, float):
                start = np.argmin(np.abs(slice_ - a.get_coord(dim)))
                updated_index_slice.append(slice(start, start+1))
            else:
                updated_index_slice.append(slice_)

        index_slice_dict = dict(zip(index_dims, updated_index_slice))
        new_slices = [slice(None) if dim not in index_dims else index_slice_dict[dim] for dim in a.dims]

        for ix, dim in enumerate(a.dims):
            a.coords[dim] = a.coords[dim][new_slices[ix]]

        a.values = a.values[tuple(new_slices)]

        return a

    def copy(self):
        return deepcopy(self)

    def merge_attrs(self, b):
        '''Merge the given dictionaries

        Args:
            b (nddata_core): attributes to merge into object
        '''

        for key in b.attrs:
            if key not in self.attrs:
                self.attrs[key] = b.attrs[key]
            else:
                if self.attrs[key] != b.attrs[key]:
                    warnings.warn('attrs in two dictionarys contain different values, leaving original value:\n{}:{}'.format(key, self.attrs[key]))

    def __len__(self):
        '''Returns total number of dims in values
        '''
        return len(self._values)

    @property
    def size(self):
        '''Return size of values
        '''
        return self._values.size

    def sort_dims(self):
        '''
        '''
        sorted_order = sorted(range(len(self.dims)), key=lambda x: self.dims[x])

#        self._dims = [self._dims[x] for x in sorted_order]
#        self._coords = [self._coords[x] for x in sorted_order]
        self.coords.reorder_index(sorted_order)
        self._values = np.moveaxis(self._values,range(len(sorted_order)),sorted_order)

    def index(self, dim):
        '''Find index of given dimension name
        '''
        if dim in self.dims:
            return self.coords.index(dim)
        else:
            raise ValueError('%s not in %s'%(dim, self.dims))

    def __truediv__(self, b):
        if isinstance(b, nddata_core):

            a, b = self.align(b)

            a.values = a.values / b.values

            # error propagation
            if a.error is not None and b.error is not None:
                error = abs(a.values) * np.sqrt((self.error/result)**2. + (b.error/result)**2.)
            elif b.error is not None:
                a.error = b.error

            return a

        elif isinstance(b, _numerical_types):
            a = self.copy()
            a.values = a.values.__truediv__(b)
            return a
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    def __rtruediv__(self, b):
        if isinstance(b, _numerical_types):
            a = self.copy()
            a.values = self.values.__rtruediv__(b)
            return a
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    def __mul__(self, b):
        if isinstance(b, nddata_core):
            a, b = self.align(b)

            a.values = a.values * b.values

            # error propagation
            if a.error is not None and b.error is not None:
                error = abs(a.values) * np.sqrt((self.error/result)**2. + (b.error/result)**2.)
            elif b.error is not None:
                a.error = b.error

            return a

        elif isinstance(b, _numerical_types):
            a = self.copy()
            a.values = a.values.__mul__(b)
            return a
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    __rmul__ = __mul__

    def __add__(self, b):
        if isinstance(b, nddata_core):
            a, b = self.align(b)

            a.values = a.values + b.values

            # error propagation
            if a.error is not None and b.error is not None:
                a.error = np.sqrt(a.error**2. + b.error**2.)
            elif b.error is not None:
                a.error = b.error

            return a

        elif isinstance(b, _numerical_types):
            a = self.copy()
            a.values = a.values.__add__(b)
            return a
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    __radd__ = __add__

    def __sub__(self, b):
        if isinstance(b, nddata_core):
            a, b = self.align(b)

            a.values = a.values - b.values

            # error propagation
            if a.error is not None and b.error is not None:
                a.error = np.sqrt(a.error**2. + b.error**2.)
            elif b.error is not None:
                a.error = b.error

            return a

        elif isinstance(b, _numerical_types):
            a = self.copy()
            a.values = a.values.__sub__(b)
            return a
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    def __rsub__(self, b):
        if isinstance(b, _numerical_types):
            a = self.copy()
            a.values = a.values.__rsub__(b)
            return a
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, b):
        if not isinstance(b, np.ndarray):
            raise TypeError('Values must be type "numpy.ndarray" not %s'%type(b))
        self._values = b

    @values.getter
    def values(self):
        return self._values

    @property
    def dims(self):
        return self.coords.dims
#        return self._dims

    @dims.setter
    def dims(self, b):
#        if not self._check_dims(b):
#            raise TypeError('dims must be list of strings')
#        self._dims = b
        self.coords.dims = b

    @property
    def coords(self):
        return self._coords

#    @coords.getter
#    def coords(self):
#        return self._coords
#
    @coords.setter
    def coords(self, b):
#        if not self._check_coords(b):
#            raise TypeError('invalid coords')
        self._coords = b

    @property
    def attrs(self):
        return self._attrs

    @attrs.setter
    def attrs(self, b):
        if isinstance(b, dict):
            self._attrs = b
        else:
            raise ValueError('attrs must be type "dict"')

    @property
    def error(self):
        return self._error

    @error.setter
    def error(self, b):
        if isinstance(b, (np.ndarray, type(None))):
            self._error = b
        else:
            raise ValueError('error must be type "numpy.ndarray"')

    def rename(self, dim, new_name):
        if dim in self.dims:
            if isinstance(new_name, str):
                self.coords.rename(dim, new_name)
#                index = self.index(dim)
#                self._dims[index] = new_name
            else:
                raise TypeError('New dimension name must be type "str" not %s'%type(new_name))
        else:
            raise ValueError('Dimension name %s is not in dims'%dim)
    

    def reorder(self, dims):
        '''
        '''

        if not self._check_dims(dims):
            raise TypeError('New dims must be list of str with no duplicates')
        for dim in dims:
            if dim not in self.dims:
                raise ValueError('no such dimension: %s'%dim)

        # Add original dims to end, remove duplicates
        dims = list(OrderedDict.fromkeys(dims + self.dims))

        new_order = [dims.index(dim) for dim in self.dims]

        self.coords.reorder(dims)
        # dims
#        self.dims = dims

        # coords
#        self.coords = [self.coords[x] for x in new_order]
#        self.coords = [self.coords[dim] for dim in dims]
#        self.coords = [self.coords[x] for x in new_order]
#        self.coords = [self.coords[new_order[x]] for x in range(len(new_order))]

        # Transpose values
        self.values = np.transpose(self.values, new_order)

    def __str__(self):
        return 'values:\n{}\ndims:\n{}\ncoords:\n{}\nattrs:\n{}'.format(self.values, self.dims, self.coords, self.attrs)

    def __repr__(self):
        return 'nddata_core(values = {}, coords = {}, dims = {}, attrs = {})'.format(repr(self.values), repr(self.dims), repr(self.coords), repr(self.attrs))

    def squeeze(self):
        '''Remove length 1 axes
        '''
        a = self.copy()
        shape = a.shape

        remove_dims = [a.dims[x] for x in range(len(shape)) if shape[x] == 1]
#
        values = np.squeeze(a.values)

        if a.error is not None:
            a.error = np.squeeze(a.error)

        attrs = a.attrs

        for dim in remove_dims:
            out = a.coords.pop(dim)
            if dim not in attrs:
                attrs[dim] = np.array(out)
            else:
                warnings.warn('Attribute lost {}:{}'.format(lost_dims[ix], lost_coords[ix]))

        return a

    def chunk(self, dim, new_dim, new_coord):
        '''
        '''

#        a = self.copy()
#        num_chunks = len(a.get_coord(dim)) / len(new_coords)
#        dims = a.dims.append(new_dim)
#        coords = a.coords.append(new_coord)
#
#        values = np.stack(np.split(a.values, num_chunks))
#        error = np.stack(np.split(a.error, num_chunks))
#
#        return a
        return NotImplemented


    def get_coord(self, dim):
        '''Return coord corresponding to given dimension name

        Args:
            dim (str): Name of dim to retrieve coordinates from

        Returns:
            numpy.ndarray: array of coordinates
            
        '''

#        return self.coords[self.index(dim)]
        return self.coords[dim]


    @property
    def shape(self):
        return self.values.shape

    @property
    def dtype(self):
        return self.values.dtype

    def sum(self, dim):
        '''Perform sum down given dimension
        '''

        a = self.copy()

        index = a.index(dim)

        a.values = a.values.sum(index)

        if a.error is not None:
            a.error = a.error.std(index)

        a.coords.pop(dim)

        return a


    def align(self, b):
        '''
        '''
        a = self.copy()
        b = b.copy()

        all_dims = list(OrderedDict.fromkeys(a.dims + b.dims))
        new_b_order = [dim for dim in all_dims if dim in b.dims]
        new_order = [b.index(dim) for dim in new_b_order]

        # create new dims where necessary
        values = a.values[tuple(slice(None) if dim in a.dims else None for dim in all_dims)]

        # re-order
        old_order = list(range(len(new_order)))
        # re-order b values so they match order of all_dims
        values_b = np.moveaxis(b.values, new_order, old_order)
        # create new dims where necessary
        values_b = values_b[tuple(slice(None) if dim in b.dims else None for dim in all_dims)]

        # Handle Error 
        if a.error is not None:
            error = a.error[tuple(slice(None) if dim in a.dims else None for dim in all_dims)]
        else: error = a.error
        if b.error is not None:
            error_b = np.moveaxis(b.values, new_order, old_order)
            error_b = error_b[tuple(slice(None) if dim in b.dims else None for dim in all_dims)]
        else: error_b = b.error

        # check coords
        for dim in all_dims:
            if (dim in a.dims) and (dim in b.dims):
                if not np.allclose(a.get_coord(dim), b.get_coord(dim)):
                    raise ValueError('Coords do not match for dim: %s'%dim)

        # merge attrs
        a.merge_attrs(b)

        # coords
        coords = a.coords.copy()
        for dim in new_b_order:
            if dim not in a.dims:
                coords.append(dim, b.coords[dim])

        a.values = values
        a.coords = coords
        a.error = error

        b.values = values_b
        b.coords = coords
        b.error = error_b

        return a, b

    def __array__(self):
        return self.values

    def sort(self, dim):
        '''
        '''
#        index = self.index(dim)

        sort_array = np.argsort(self.coords[dim])

#        self.coords[index] = self.coords[index][sort_array]
        self.coords[dim] = self.coords[dim][sort_array]

        new_order = tuple([slice(None) if dim != this_dim else sort_array for this_dim in self.dims])

        self.values = self.values[new_order]

    def is_sorted(self, dim):
        '''
        '''
        return np.all(self.coords[dim][:-1] <= self.coords[dim][1:])

    @property
    def real(self):
        a = self.copy()
        a.values = np.real(a.values)
        return a

    @property
    def imag(self):
        a = self.copy()
        a.values = np.imag(a.values)
        return a

    @property
    def abs(self):
        a = self.copy()
        a.values = np.abs(a.values)
        return a

    def concatenate(self, b, dim):
        '''
        '''

        if not dim in b.dims:
            raise ValueError('dim does not exist')
        if not dim in self.dims:
            raise ValueError('dim does not exist')

        index = self.dims.index(dim)

        b.reorder(self.dims)

        self.values = np.concatenate((self.values, b.values), axis = index)
        
        self.coords[dim] = np.concatenate((np.array(self.coords[dim]).reshape(-1),np.array(b.coords[dim]).reshape(-1)))

    def new_dim(self, dim, coord):
        '''
        '''

        self.coords.append(dim, np.r_[coord])
        self.values = np.expand_dims(self.values, -1)


    def argmax(self, dim):
        '''Return argmax for given dim
        '''
        a = self.copy()
        index = a.dims.index(dim)

        a.values = a.coords[dim][np.argmax(a.values, axis = index)]
        a.coords.pop(dim)

        return a

    def argmin(self, dim):
        '''Return argmin for given dim
        '''
        a = self.copy()
        index = a.dims.index(dim)

        a.values = a.coords[dim][np.argmin(a.values, axis = index)]
        a.coords.pop(dim)

        return a

    @property
    def ndim(self):
        return self.values.ndim

if __name__ == '__main__':
#    a = np.array(range(9)).reshape(3,3)
#    b = ['x','y']
#    c = [np.r_[[1,2,3]],np.r_[[4,5,6]]]
#    data = nddata_core(np.array(range(27)).reshape(3,3,3), [np.r_[1,2,3],np.r_[4,5,6], np.r_[7,8,9]], ['x','y','z'])
#    data = nddata_core(np.array(range(27)).reshape(3,3,3), ['x','y','z'], [np.r_[1,2,3],np.r_[4,5,6], np.r_[7,8,9]])

    x = np.r_[0:10:.1]
    y = np.r_[0:20]
    z = np.r_[0:15]
    random_order = np.argsort(np.random.randn(len(x)))
    x = x[random_order]

#    q = np.r_[0:5]
#    r = np.r_[0:17]
#    p = np.r_[0:13]
#    data = nddata_core(np.array(range(len(x)*len(y)*len(p))).reshape(len(x),len(y),len(p)), ['x','y','p'], [x, y, p])
    data = nddata_core(np.array(range(len(x)*len(y)*len(z))).reshape(len(x),len(y),len(z)), ['x','y','z'], [x, y, z])
    data.reorder(['x','y','z'])

    d = nddata_core(x, ['x'], [x])
    d2 = data['x',0]
    d3 = d2.squeeze()
