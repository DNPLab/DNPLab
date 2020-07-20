from __future__ import division
import operator
import numpy as np
import warnings
from copy import deepcopy

class nddata:
    '''nddata class
    '''

    def __init__(self, values = np.r_[[]], dims = [], coords = [], attrs = {}, error = None, **kwargs):

        if isinstance(values, np.ndarray):
            self._values = values
        else:
            raise TypeError('values must be type "numpy.ndarray" not %s'%str(type(values)))

        if self._check_coords(coords):
            self._coords = coords
        else:
            raise TypeError('Invalid coords. Must be list of numpy.ndarray')

        if self._check_dims(dims):
            self._dims = dims
        else:
            raise TypeError('Invalid dims. Must be list of str with no duplicates')

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
            raise ValueError('Dimensions not consistent')

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
            coords_check = list(self._values.shape) == [len(self._coords[x]) for x in range(len(self._coords))]

        dims_check = len(self.values.shape) == len(self.dims)

        return coords_check and dims_check

    def _attrs_valid():
        '''Verify attrs are valid
        '''

        for key in self._attrs:
            if not isinstance(self._attrs, (list, np.ndarray, int, float, str)):
                return False

        return True

    def copy(self):
        return deepcopy(self)

    def merge_attrs(self, a, b):
        '''Merge the given dictionaries

        Args:
            a (dict): Starting Dictionary
            b (dict): Dictionary to merge into attributes
        '''
        if not isinstance(b, dict):
            raise ValueError('Must be dict')

        for key in b:
            if key not in a:
                a[key] = b[key]
            else:
                if a[key] != b[key]:
                    warnings.warn('attrs in two dictionarys contain different values, leaving original value:\n{}:{}'.format(key, self._attrs[key]))
        return a

    def __len__(self):
        '''Returns total number of dims in values
        '''
        return len(self._values)

    def size(self):
        '''Return size of values
        '''
        return self._values.size

    def sort(self):
        '''
        '''
        sorted_order = sorted(range(len(self._dims)), key=lambda x: self._dims[x])

        self._dims = [self._dims[x] for x in sorted_order]
        self._coords = [self._coords[x] for x in sorted_order]
        self._values = np.moveaxis(self._values,range(len(sorted_order)),sorted_order)


    def index(self, dim):
        '''Find index of given dimension name
        '''
        if dim in self.dims:
            return self.dims.index(dim)
        else:
            raise ValueError('%s not in %s'%(dim, self.dims))

    def __truediv__(self, b):
        if isinstance(b, nddata):
            # reorder 
            old_order = b.dims
            b.reorder(self.dims)

            # check coords
            for ix in range(len(self.coords)):
                if not np.allclose(self.coords[ix],b.coords[ix]):
                    raise ValueError('Coords do not match')

            # add result
            result = self.values.__truediv__(b.values)

            # merge attrs
            attrs = self.merge_attrs(self.attrs, b.attrs)

            # error propagation
            if self._check_error(b.error):
                if self.error is not None:
                    error = abs(result) * np.sqrt((self.error/result)**2. + (b.error/result)**2.)
                else:
                    error = None
            else:
                error = self.error

            # return b to original order
            b.reorder(old_order)

            return nddata(result, self.dims, self.coords, attrs, error)

        elif isinstance(b, (np.ndarray, int, float)):
            result = self.values.__truediv__(b)
            return nddata(result, self.dims, self.coords, self.attrs, self.error)
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    def __rtruediv__(self, b):
        if isinstance(b, (np.ndarray, int, float)):
            result = self.values.__rtruediv__(b)
            return nddata(result, self.dims, self.coords, self.attrs, self.error)
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    def __mul__(self, b):
        if isinstance(b, nddata):
            # reorder 
            old_order = b.dims
            b.reorder(self.dims)

            # check coords
            for ix in range(len(self.coords)):
                if not np.allclose(self.coords[ix],b.coords[ix]):
                    raise ValueError('Coords do not match')

            # add result
            result = self.values.__mul__(b.values)

            # merge attrs
            attrs = self.merge_attrs(self.attrs, b.attrs)

            # error propagation
            if self._check_error(b.error):
                if self.error is not None:
                    error = abs(result) * np.sqrt((self.error/result)**2. + (b.error/result)**2.)
                else:
                    error = None
            else:
                error = self.error

            # return b to original order
            b.reorder(old_order)

            return nddata(result, self.dims, self.coords, attrs, error)

        elif isinstance(b, (np.ndarray, int, float)):
            result = self.values.__mul__(b)
            return nddata(result, self.dims, self.coords, self.attrs, self.error)
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    __rmul__ = __mul__

    def __add__(self, b):
        if isinstance(b, nddata):
            # reorder 
            old_order = b.dims
            b.reorder(self.dims)

            # check coords
            for ix in range(len(self.coords)):
                if not np.allclose(self.coords[ix],b.coords[ix]):
                    raise ValueError('Coords do not match')

            # add result
            result = self.values + b.values

            # merge attrs
            attrs = self.merge_attrs(self.attrs, b.attrs)

            # error propagation
            if self._check_error(b.error):
                if self._error is not None:
                    error = np.sqrt(self.error**2. + b.error**2.)
                else:
                    error = None
            else:
                error = self.error

            # return b to original order
            b.reorder(old_order)

            return nddata(result, self.dims, self.coords, attrs, error)

        elif isinstance(b, (np.ndarray, int, float)):
            result = self.values + b
            return nddata(result, self.dims, self.coords, self.attrs, self.error)
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    __radd__ = __add__

    def __sub__(self, b):
        if isinstance(b, nddata):
            # reorder 
            old_order = b.dims
            b.reorder(self.dims)

            # check coords
            for ix in range(len(self.coords)):
                if not np.allclose(self.coords[ix],b.coords[ix]):
                    raise ValueError('Coords do not match')

            # add result
            result = self.values.__sub__(b.values)

            # merge attrs
            attrs = self.merge_attrs(self.attrs, b.attrs)

            # error propagation
            if self._check_error(b.error):
                if self._error is not None:
                    error = np.sqrt(self.error**2. + b.error**2.)
                else:
                    error = None
            else:
                error = self.error

            # return b to original order
            b.reorder(old_order)

            return nddata(result, self.dims, self.coords, attrs, error)

        elif isinstance(b, (np.ndarray, int, float)):
            result = self.values.__sub__(b)
            return nddata(result, self.dims, self.coords, self.attrs, self.error)
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    def __rsub__(self, b):
        if isinstance(b, (np.ndarray, int, float)):
            result = self.values.__rsub__(b)
            return nddata(result, self.dims, self.coords, self.attrs, self.error)
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
        return self._dims

    @dims.setter
    def dims(self, b):
        if not self._check_dims(b):
            raise TypeError('dims must be list of strings')
        self._dims = b

    @property
    def coords(self):
        return self._coords

    @coords.getter
    def coords(self):
        return self._coords

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
        if isinstance(b, np.ndarray):
            self._error = b
        else:
            raise ValueError('error must be type "numpy.ndarray"')

    def rename(self, dim, new_name):
        if dim in self._dims:
            if isinstance(new_name, str):
                index = self.index(dim)
                self._dims[index] = new_name
            else:
                raise TypeError('New dimension name must be type "str" not %s'%type(new_name))
        else:
            raise ValueError('Dimension name %s is not in dims'%dim)
    
    def align(self, a, b):
        '''Align two nddata objects

        Args:
            a (nddata): First nddata object, these dimensions are given first
            b (nddata): Second nddata object
        '''
        if not isinstance(a, nddata):
            raise ValueError('Must be type nddata not %s'%str(type(a)))
        if not isinstance(b, nddata):
            raise ValueError('Must be type nddata not %s'%str(type(b)))

        all_dims = list(set(a.dims + b.dims))

        # Make sure coords match when dims match
        for dim in all_dims:
            if (dim in a.dims) and (dim in b.dims):
                if not np.allclose(a.get_coord[dim], b.get_coord[dim]):
                    raise ValueError('Axes')

        new_a = a.copy()
        new_b = b.copy()

        new_a = new_a.reorder(all_dims)
        new_b = new_b.reorder(all_dims)

        return new_a, new_b

    def operate(self, a, b, operation):
        '''
        '''

        return operation(a, b)


    def reorder(self, new_dims):
        '''
        '''

        if not self._check_dims(new_dims):
            raise TypeError('New dims must be list of str with no duplicates')
        if not len(self._dims) == len(new_dims):
            raise ValueError('Number of dims do not match')

        new_order = [new_dims.index(dim) for dim in self._dims]
        print(new_order)

        self._dims = new_dims#[self._dims[x] for x in new_order]
        self._coords = [self._coords[x] for x in new_order]
        self._values = np.moveaxis(self._values,range(len(new_order)),new_order)

    def __str__(self):
        return 'values:\n{}\ndims:\n{}\ncoords:\n{}\nattrs:\n{}'.format(self._values, self._dims, self._coords, self._attrs)

    def __repr__(self):
        return 'nddata(values = {}, coords = {}, dims = {}, attrs = {})'.format(repr(self._values), repr(self._dims), repr(self._coords), repr(self._attrs))

    def squeeze(self):
        '''Remove length 1 axes
        '''
        shape = self.shape

        dims = [self.dims[x] for x in range(len(shape)) if shape[x] != 1]
        coords = [self.coords[x] for x in range(len(shape)) if shape[x] != 1]

        lost_dims = [self.dims[x] for x in range(len(shape)) if shape[x] == 1]
        lost_coords = [self.coords[x] for x in range(len(shape)) if shape[x] == 1]

        values = np.squeeze(self.values)
        error = np.squeeze(self.error)

        attrs = self.attrs

        for ix in range(len(lost_dims)):
            if lost_dims[ix] not in attrs:
                attrs[lost_dims[ix]] = lost_coords[ix]
            else:
                warnings.warn('Attribute lost {}:{}'.format(lost_dims[ix],lost_coords[ix]))

        return nddata(values, dims, coords, attrs, error)

    def chunk(self, dim, new_dim, new_coord):
        '''
        '''

        num_chunks = len(self.get_coord(dim)) / len(new_coords)
        dims = self.dims.append(new_dim)
        coords = self.coords.append(new_coord)

        values = np.stack(np.split(self.values, num_chunks))
        error = np.stack(np.split(self.error, num_chunks))

        return nddata(values, dims, coords, attrs, error)


    def get_coord(self, dim):
        '''Return coord corresponding to given dimension name

        Args:
            dim (str): Name of dim to retrieve coordinates from

        Returns:
            numpy.ndarray: array of coordinates
            
        '''

        return self.coords[self.index(dim)]


    @property
    def shape(self):
        return self.values.shape

if __name__ == '__main__':
#    a = np.array(range(9)).reshape(3,3)
#    b = ['x','y']
#    c = [np.r_[[1,2,3]],np.r_[[4,5,6]]]
#    data = nddata(np.array(range(27)).reshape(3,3,3), [np.r_[1,2,3],np.r_[4,5,6], np.r_[7,8,9]], ['x','y','z'])
#    data = nddata(np.array(range(27)).reshape(3,3,3), ['x','y','z'], [np.r_[1,2,3],np.r_[4,5,6], np.r_[7,8,9]])

    x = np.r_[1:3]
    y = np.r_[1:4]
    z = np.r_[1]
    data = nddata(np.array(range(len(x)*len(y)*len(z))).reshape(len(x),len(y),len(z)), ['x','y','z'], [x, y, z])
    print(data)
