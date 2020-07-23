from __future__ import division
import operator
import numpy as np
import warnings
from copy import deepcopy
from collections import OrderedDict

class nddata_core(object):
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
            warnings.warn('Dimensions not consistent')
#            raise ValueError('Dimensions not consistent')

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

    @property
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
        if isinstance(b, nddata_core):
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

            return nddata_core(result, self.dims, self.coords, attrs, error)

        elif isinstance(b, (np.ndarray, int, float)):
            result = self.values.__truediv__(b)
            return nddata_core(result, self.dims, self.coords, self.attrs, self.error)
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    def __rtruediv__(self, b):
        if isinstance(b, (np.ndarray, int, float)):
            result = self.values.__rtruediv__(b)
            return nddata_core(result, self.dims, self.coords, self.attrs, self.error)
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    def __mul__(self, b):
        if isinstance(b, nddata_core):
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

            return nddata_core(result, self.dims, self.coords, attrs, error)

        elif isinstance(b, (np.ndarray, int, float)):
            result = self.values.__mul__(b)
            return nddata_core(result, self.dims, self.coords, self.attrs, self.error)
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    __rmul__ = __mul__

    def __add__(self, b):
        if isinstance(b, nddata_core):

            all_dims = list(OrderedDict.fromkeys(self.dims + b.dims))
            new_b_order = [dim for dim in all_dims if dim in b.dims]
            new_order = [b.index(dim) for dim in new_b_order]

            # create new dims where necessary
            values = self.values[tuple(slice(None) if dim in self.dims else None for dim in all_dims)]

            # re-order
            old_order = list(range(len(new_order)))
            # re-order b values so they match order of all_dims
            values_b = np.moveaxis(b.values, new_order, old_order)
            # create new dims where necessary
            values_b = values_b[tuple(slice(None) if dim in b.dims else None for dim in all_dims)]

            # Handle Error 
            if self.error is not None:
                error = self.error[tuple(slice(None) if dim in self.dims else None for dim in all_dims)]
            if b.error is not None:
                error_b = np.moveaxis(b.values, new_order, old_order)
                error_b = error_b[tuple(slice(None) if dim in b.dims else None for dim in all_dims)]

            # check coords
            for dim in all_dims:
                if (dim in self.dims) and (dim in b.dims):
                    if not np.allclose(self.get_coord(dim), b.get_coord(dim)):
                        raise ValueError('Coords do not match for dim: %s'%dim)

            # add result
            result = values + values_b

            # merge attrs
            attrs = self.merge_attrs(self.attrs, b.attrs)

            # error propagation
            if self._check_error(b.error):
                if self._error is not None:
                    error = np.sqrt(error**2. + error_b**2.)
                else:
                    error = None
            else:
                error = self.error

            coords = list(self.coords)
            coords += [b.coords[ix] for ix in new_order if b.dims[ix] not in self.dims]

            return nddata_core(result, all_dims, coords, attrs, error)

        elif isinstance(b, (np.ndarray, int, float)):
            result = self.values + b
            return nddata_core(result, self.dims, self.coords, self.attrs, self.error)
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    __radd__ = __add__

    def __sub__(self, b):
        if isinstance(b, nddata_core):
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

            return nddata_core(result, self.dims, self.coords, attrs, error)

        elif isinstance(b, (np.ndarray, int, float)):
            result = self.values.__sub__(b)
            return nddata_core(result, self.dims, self.coords, self.attrs, self.error)
        else:
            raise TypeError('Cannot add type: {}'.format(type(b)))

    def __rsub__(self, b):
        if isinstance(b, (np.ndarray, int, float)):
            result = self.values.__rsub__(b)
            return nddata_core(result, self.dims, self.coords, self.attrs, self.error)
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
    
    def _align(self, b):
        '''Align two nddata objects

        Args:
            a (nddata): First nddata object, these dimensions are given first
            b (nddata): Second nddata object
        '''
        if not isinstance(self, nddata_core):
            raise ValueError('Must be type nddata not %s'%str(type(self)))
        if not isinstance(b, nddata_core):
            raise ValueError('Must be type nddata not %s'%str(type(b)))

#        all_dims = list(set(self.dims + b.dims)) # does not keep order
        # slow method
        all_dims = list(OrderedDict.fromkeys(self.dims+b.dims))
        # fast method here: https://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-whilst-preserving-order


        # Make sure coords match when dims match
        for dim in all_dims:
            if (dim in self.dims) and (dim in b.dims):
                if not np.allclose(self.get_coord(dim), b.get_coord(dim)):
                    raise ValueError('Coords do not match for dim %s'%dim)

        dims = all_dims
        coords = self.coords
        coords += [b.get_coord(dim) for dim in b.dims if dim not in self.dims]

        b_dims = all_dims
        b_coords = coords

#        print(b.values[(slice(None) if dim in b.dims else None for dim in all_dims)])

#        print(all_dims)
#        print(dims)
#        print(coords)

        new_shape = tuple(self.get_coord(dim).shape[0] if dim in self.dims else 1 for dim in all_dims)

        # Must re-order not reshape

#        b_new_shape = tuple(b.get_coord(dim).shape[0] if dim in b.dims else 1 for dim in all_dims)
#        print(new_shape)

        values = self.values.reshape(new_shape)
        if self.error is not None:
            error = self.error.reshape(new_shape)
#        add_axis = (ix for range()

        return nddata_core(values, dims, coords)

    def operate(self, a, b, operation):
        '''
        '''

        return operation(a, b)


    def reorder(self, new_dims, force_new_dims = False):
        '''
        '''

        if not force_new_dims:
            if not self._check_dims(new_dims):
                raise TypeError('New dims must be list of str with no duplicates')
            if not len(self._dims) == len(new_dims):
                raise ValueError('Number of dims do not match')

        # old method, no force new dims
#        new_order = [new_dims.index(dim) for dim in self._dims]

        new_order = [new_dims.index(dim) for dim in self.dims if new_dims in self.dims]

        # reshape 
        print(new_order)

        self._dims = new_dims#[self._dims[x] for x in new_order]

        # not forced
#        self._coords = [self.coords[x] for x in new_order]
        #forced
        self._coords = [self.coords.index(new_dims[x]) if new_dims[x] in self.dims else np.r_[[]] for x in range(len(new_dims))]

#        self._values

        # reorder axes if they exist
        self._values = np.moveaxis(self._values,range(len(new_order)),new_order)

        # add new dimensions if forced values
        self._values[(slice(None) if dim in b.dims else None for dim in all_dims)]

    def __str__(self):
        return 'values:\n{}\ndims:\n{}\ncoords:\n{}\nattrs:\n{}'.format(self._values, self._dims, self._coords, self._attrs)

    def __repr__(self):
        return 'nddata_core(values = {}, coords = {}, dims = {}, attrs = {})'.format(repr(self._values), repr(self._dims), repr(self._coords), repr(self._attrs))

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

        return nddata_core(values, dims, coords, attrs, error)

    def chunk(self, dim, new_dim, new_coord):
        '''
        '''

        num_chunks = len(self.get_coord(dim)) / len(new_coords)
        dims = self.dims.append(new_dim)
        coords = self.coords.append(new_coord)

        values = np.stack(np.split(self.values, num_chunks))
        error = np.stack(np.split(self.error, num_chunks))

        return nddata_core(values, dims, coords, attrs, error)


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

    def sum(self, dim):
        '''Perform sum down given dimension
        '''

        index = self.index(dim)

        values = self.values.sum(index)

        if self.error is not None:
            error = self.error.std(index)
        else:
            error = self.error

        dims = list(self.dims)
        dims.pop(index)
        coords = list(self.coords)
        coords.pop(index)

        attrs = self.attrs

        return nddata_core(values, dims, coords, attrs, error)


    def einsum(self, b):
        '''
        '''

        all_dims = list(OrderedDict.fromkeys(self.dims+b.dims))

#        ein_dims = [ix for ix in range(len(self.dims)) if self.dims[ix] in all_dims]
#        ein_dims_b = [ix for ix in range(len(b.dims)) if b.dims[ix] in all_dims]

        ein_dims = [all_dims.index(dim) for dim in all_dims if dim in self.dims]
        ein_dims_b = [all_dims.index(dim) for dim in all_dims if dim in b.dims]
        print(ein_dims)
        print(ein_dims_b)

        values = np.einsum(self.values, ein_dims, b.values, ein_dims_b)

        dims = all_dims
        coords = self.coords
        coords += [b.coords[ix] for ix in range(len(b.coords)) if b.dims[ix] not in self.dims]

        attrs = self.attrs
        error = self.error

        return nddata_core(values, dims, coords, attrs, error)

#    def new_dims(self, b, all_dims):
    def new_dims(self, b):
        '''
        '''


        all_dims = list(OrderedDict.fromkeys(self.dims + b.dims))

        new_shape = tuple(slice(None) if dim in self.dims else None for dim in all_dims)

        values = self.values[new_shape]
        new_shape_b = tuple(slice(None) if dim in b.dims else None for dim in all_dims)
        values_b = b.values[new_shape_b]

        return values

    def _operate_on_values(self, b, o):
        '''
        '''

        all_dims = list(OrderedDict.fromkeys(self.dims + b.dims))
        new_b_order = [dim for dim in all_dims if dim in b.dims]
        new_order = [b.index(dim) for dim in new_b_order]

        # create new dims where necessary
        values = self.values[tuple(slice(None) if dim in self.dims else None for dim in all_dims)]

        # re-order
        old_order = list(range(len(new_order)))
        # re-order b values so they match order of all_dims
        values_b = np.moveaxis(b.values, new_order, old_order)
        # create new dims where necessary
        values_b = values_b[tuple(slice(None) if dim in b.dims else None for dim in all_dims)]

        result = o(values, values_b)

        return result

    def _prepare_operation(self, b):
        '''

        Returns:
            values
            values_b
            error
            error_b
            dims
            coords
            attrs
        '''
        all_dims = list(OrderedDict.fromkeys(self.dims + b.dims))
        new_b_order = [dim for dim in all_dims if dim in b.dims]
        new_order = [b.index(dim) for dim in new_b_order]

        # create new dims where necessary
        values = self.values[tuple(slice(None) if dim in self.dims else None for dim in all_dims)]

        # re-order
        old_order = list(range(len(new_order)))
        # re-order b values so they match order of all_dims
        values_b = np.moveaxis(b.values, new_order, old_order)
        # create new dims where necessary
        values_b = values_b[tuple(slice(None) if dim in b.dims else None for dim in all_dims)]

        result = o(values, values_b)

        return values, values_b, error, error_b, dims, coords, attrs


if __name__ == '__main__':
#    a = np.array(range(9)).reshape(3,3)
#    b = ['x','y']
#    c = [np.r_[[1,2,3]],np.r_[[4,5,6]]]
#    data = nddata_core(np.array(range(27)).reshape(3,3,3), [np.r_[1,2,3],np.r_[4,5,6], np.r_[7,8,9]], ['x','y','z'])
#    data = nddata_core(np.array(range(27)).reshape(3,3,3), ['x','y','z'], [np.r_[1,2,3],np.r_[4,5,6], np.r_[7,8,9]])

    x = np.r_[0:10]
    y = np.r_[0:20]
    z = np.r_[0:15]
    q = np.r_[0:5]
    r = np.r_[0:17]
    p = np.r_[0:13]
    data = nddata_core(np.array(range(len(x)*len(y)*len(p))).reshape(len(x),len(y),len(p)), ['x','y','p'], [x, y, p])

    data2 = nddata_core(np.array(range(len(x)*len(y)*len(q))).reshape(len(x),len(y),len(q)), ['x','y','q'], [x, y, q])

    data3 = nddata_core(np.array(range(len(r)*len(p)*len(q))).reshape(len(r),len(p),len(q)), ['r','p','q'], [r, p, q])

    d = data + data
    d = data + data2
    d = data + data3

    d = data2 + data2
    print('-'*50)
    d = data2 + data3

    print(data._operate_on_values(data3,operator.__add__))
    print(data._operate_on_values(data3,operator.__add__).shape)




#    d = data + data3
#    print(data)

#    data._align(data2)
#    data.einsum(data2)
#    data.new_dims(['blue','foo'])
#    data.new_dims(data2)
#    print(data + data2)

