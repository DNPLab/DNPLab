from __future__ import division
import numpy as np
import operator
import functools
from collections import OrderedDict
from copy import deepcopy


#allowed_domains = ['FT','IFT','LT','ILT','Wavelet','IWavelet']
allowed_domains = ['FT','IFT']


class nddata_coord(object):
    '''
    '''

    def __init__(self, dim, *args, step_type = 'linear', endpoint = False,**kwargs):
        '''
        Args:
            dim (str): name of dimension
        '''

        self.dim = dim
        if len(args) == 0:
            self.__start = None
            self.__stop = None
            self.__step = None

        if len(args) == 1:
            coord = args[0]
            if isinstance(coord, np.ndarray):
                self._axis = coord
                try:
                    self.reduce()
                    self._type = 'linear'
                except:
                    self._type = 'ndarray'
            elif isinstance(coord, slice):
                self.__start = coord.start
                self.__stop = coord.stop
                self.__step = coord.step
            elif isinstance(coord, (int,float)):
                self.__start = None
                self.__stop = coord
                self.__step = None
            else:
                raise ValueError('coord not understood')
        if len(args) == 2:
            start = args[0]
            stop = args[1]
            if isinstance(start, (int, float)) and isinstance(stop, (int,float)):
                self.__start = start
                self.__stop = stop
                self.__step = None
            else:
                raise TypeError('coord not understood')
        if len(args) == 3:
            start = args[0]
            stop = args[1]
            step = args[2]
            if isinstance(start, (int,float)) and isinstance(stop, (int, float)) and isinstance(step, (int, float)):
                self.__start = start
                self.__stop = stop
                self.__step = step
            else:
                raise TypeError('coord not understood')

        # set up other variables
        if 'type' in kwargs:
            self._type = kwargs['type']
        else:
            self._type = 'linear'
    
    def transform(self, new_domain, shift = False):
        '''
        '''

        if new_domain not in allowed_domains:
            raise ValueError('domain not supported')

        if new_domain == 'FT' or new_domain == 'IFT':
            step = 1. / (self.size * self.step) # Does this reduce performance? size calculates array

            if shift == True:
                start = -0.5/self.step
                stop = 0.5/self.step
            else:
                start = 0
                stop = 1./self.step

        if hasattr(self, '_array'):
            del self._array

        return nddata_coord(self.dim, slice(start, stop, step))
        
    def reduce(self):
        '''
        '''

        if hasattr(self, '_array'):
            start = self._array[0]
            step = self._array[1] - self._array[0]
            stop = self._array[-1] + step

            if np.allclose(self._array, np.r_[slice(start,stop,step)]):
                self.start = start
                self.stop = stop
                self.step = step
            else:
                raise ValueError('Array must have evenly spaced values to be reduced')
        else:
            raise ValueError('No array to reduce')

    @property
    def dim(self):
        return self.__dim

    @dim.setter
    def dim(self, dim):
        '''
        '''

        if isinstance(dim, str):
            self.__dim = dim
        else:
            raise TypeError('dim must be type str not %s'%str(type(dim)))

    @property
    def start(self):
        '''
        '''
        return self.__start

    @start.setter
    def start(self, b):
        '''
        '''
        if isinstance(b, (int, float)):
            self.__start = b
        if hasattr(self, '_array'):
            del self._array

    @property
    def stop(self):
        '''
        '''
        return self.__stop

    @stop.setter
    def stop(self, b):
        '''
        '''
        if isinstance(b, (int, float)):
            self.__stop = b
        if hasattr(self, '_array'):
            del self._array

    @property
    def step(self):
        '''
        '''
        return self.__step

    @step.setter
    def step(self, b):
        '''
        '''
        if isinstance(b, (int, float)):
            self.__step = b
        if hasattr(self, '_array'):
            del self._array

    @property
    def array(self):
        '''Return axes as numpy array
        '''

        # Fast method
#        try:
#            return self._array
#        except:
#            self._array = np.r_[slice(self.start, self.stop, self.step)]
#            return self._array
        # 2nd Attempt at fast method
        if hasattr(self, '_array'):
            return self._array
        else:
            self._array = np.r_[slice(self.start, self.stop, self.step)]
            return self._array

#        # Slow method, 56 times slower
#        return np.r_[slice(self.start, self.stop, self.step)]

    @array.setter
    def array(self, b):
        '''
        '''
        if isinstance(b, np.ndarray):
            self._array = b
        else:
            raise TypeError('array type must be numpy.ndarray')
           

    @array.deleter
    def array(self):
        del self._array

    def slice(self, *args):
        '''
        '''
        return self.array[slice(*args)]

    def __getitem__(self, x):
        '''
        '''
        return self.array[x] # Faster by 10% if array is stored in object 
#        return self.start + self.step * x

    def _del_array(self):
        '''
        '''
        del self._array


    @property
    def size(self):
        '''
        '''
        if hasattr(self, '_array'):
            return self.array.size # 3 times faster when array is stored in object
        else:
            return max(int((self.stop - self.start) / self.step),0) # Take max for case when direction is reversed

    def __repr__(self):
        return 'nddata_coord(\'{}\', {})'.format(self.dim,self.array)

    def __str__(self):
        return '\'{}\':{}'.format(self.dim,str(self.array))
    
    def __add__(self, b):
        start = self.start + b
        stop = self.stop + b
        step = self.step

        if hasattr(self, '_array'):
            del self.array

        return nddata_coord(self.dim, slice(start, stop, step))

    __radd__ = __add__

    def __sub__(self, b):
        start = self.start - b
        stop = self.stop - b
        step = self.step

        if hasattr(self, '_array'):
            del self.array

        return nddata_coord(self.dim, slice(start, stop, step))

    def __rsub__(self, b):
        # switch start and stop
        start = b - self.stop
        stop = b - self.start 
        step = self.step

        if hasattr(self, '_array'):
            del self.array

        return nddata_coord(self.dim, slice(start, stop, step))

    def __mul__(self, b):
        start = b * self.start
        stop = b * self.stop
        step = b * self.step

        if hasattr(self, '_array'):
            del self.array

        return nddata_coord(self.dim, slice(start, stop, step))

    __rmul__ = __mul__

    def __truediv__(self, b):
        start = self.start / b
        stop = self.stop / b
        step = self.step / b

        if hasattr(self, '_array'):
            del self.array

        return nddata_coord(self.dim, slice(start, stop, step))

    def __rtruediv__(self, b):
        start = b / self.start
        stop = b / self.stop
        step = b / self.step

        if hasattr(self, '_array'):
            del self.array

        return nddata_coord(self.dim, slice(start, stop, step))

    def __array__(self):
        return self.array

    def __len__(self):
        return len(self.array)

    def __matmul__(self, b):
        return self.start + self.step * b

    @property
    def shape(self):
        return self.array.shape

#class nddata_coord_collection(MutableMapping):
class nddata_coord_collection(object):
    def __init__(self, dims, coords):

        if not isinstance(dims, list):
            raise ValueError('dims must be a list')
        if not isinstance(coords, list):
            raise ValueError('coords must be a list')

        if self._check_dims:
            self._dims = dims
        else:
            raise TypeError('dims must be list of str')

        if self._check_coords:
            self._coords = coords
        else:
            raise TypeError('coords must be list of 1d numpy arrays')

#        for arg in args:
#            if isinstance(arg, nddata_coord):
#                self.__dims.append(arg.dim)
#                self.__coords.append(arg)
#            elif isinstance(arg, np.ndarray):
#                raise ValueError('No dim given for numpy ndarray')
#
#        for kwarg in kwargs: #key, value pairs only in order for python >=3.6
#            dim = kwarg
#            coord = kwarg[dim]
#            if isinstance(coord, np.ndarray):
#                self.__dims.append(kwarg)
#                self.__coords.append(nddata_coord(dim,coord))
#            else:
#                raise TypeError('kwarg not understood')

    def _check_dims(self, dims):
        '''Verify dims is a list of str

        Args:
            dims: Object to test for valid dims

        Returns:
            bool: True if valid dims, False otherwise
        '''
        # Check if dims is list
        if not isinstance(dims, list):
            return False
        
        # Check if all members are str
        for dim in dims:
            if not isinstance(dim, str):
                return False

        # Check for duplicates
        if len(dims) != len(set(dims)):
            return False

        return True

    def _check_coords(self, coords):
        '''Check if valid coords

        Args:
            coords: Object to test

        Returns:
            bool: True if list of 1d numpy arrays, False otherwise
        '''

        # Verify coords is a list
        if not isinstance(coords, list):
            return False

        # Verify each member is 1d numpy array (and not empty)
        for coord in coords:
            if not isinstance(coord, (nddata_coord, np.ndarray)):
                return False
            if (not (len(coord.shape) == 1)) or (coord.size == 0):
                return False

        return True

    def _self_consistent(self):
        return NotImplemented

    def index(self, dim):
        return self._dims.index(dim)

    def __getitem__(self, dim):
        if isinstance(dim, str):
            return self.coords[self.index(dim)]
        elif isinstance(dim, int):
            return self.coords[dim]
        else:
            raise TypeError('dim must be type str or int not: %s'%str(type(dim)))

    def __setitem__(self, dim, coord):
        if not isinstance(dim, str):
            raise TypeError('dim must be type str not %s'%str(type(dim)))
        if not isinstance(coord, (nddata_coord, np.ndarray)):
            raise TypeError('argument must be type nddata_coord or numpy ndarray not %s'%str(type(coord)))

#        if isinstance(coord, np.ndarray):
#            coord = nddata_coord(dim, coord)

        # if dim already in dims, overwrite
        if dim in self.dims:
            index = self.index(dim)
            self._coords[index] = coord
        else:
            self._dims.append(dim)
            self._coords.append(coord)

    def __delitem__(self, dim):
        index = self.index(dim)
        del self._dims[index]
        del self._coords[index]

    @property
    def dims(self):
        return self._dims

    @dims.setter
    def dims(self, dims):
        if self._check_dims(dims):
            self._dims = dims
        else:
            raise TypeError('Invalid dims. Cannot set dims to {}'.format(dims))

    @property
    def coords(self):
        return np.array(self._coords)

    @coords.setter
    def coords(self, coords):
        if self._check_coords(coords):
            self._coords = coords
        else:
            raise TypeError('Invalid coords. Cannot set coords to {}'.format(coords))

    def __repr__(self):
        return 'nddata_coord_collection({})'.format(self.coords)

    def __str__(self):
        return 'dims:\n{}\ncoords:\n{}'.format(self.dims, self.coords)

    def __iter__(self):
        return iter(self.coords)

    def __len__(self):
        return len(self.coords)

    @property
    def shape(self):
        return tuple(k.size for k in self.coords)

    @property
    def size(self):
        return functools.reduce(operator.mul, [len(k) for k in self.coords], 1)

    def reorder(self, dims):

        if not self._check_dims(dims):
            raise TypeError('New dims must be list of str with no duplicates')
        for dim in dims:
            if dim not in self.dims:
                raise ValueError('no such dimension: %s'%dim)

        # Add original dims to end, remove duplicates
        dims = list(OrderedDict.fromkeys(dims + self.dims))

        # New indices for dims
        new_order = [dims.index(dim) for dim in self.dims]

        self.dims = dims

        self.coords = [self.coords[x] for x in new_order]

    def pop(self, dim):
        index = self.index(dim)
        
        out = self._coords.pop(index)

        self.dims.pop(index)

        return out

    def copy(self):
        '''
        '''
        return deepcopy(self)
    def reorder_index(self, new_order):
        '''Reorder based on index
        '''

        self._coords = [self._coords[x] for x in new_order]
        self._dims = [self._dims[x] for x in new_order]

    def rename(self, dim, new_dim):
        '''
        '''

        if isinstance(self[dim], nddata_coord):
            self[dim].dim = new_dim
            self.dims[self.index(dim)] = new_dim
        else:
            self.dims[self.index(dim)] = new_dim

    def append(self, dim, coord):
        '''Append to coords
        '''
        if not isinstance(dim, str):
            raise TypeError('dim must be type str not %s'%str(type(dim)))

        if not isinstance(coord, (complex, float, int, np.ndarray, nddata_coord)):
            raise TypeError('coord must be type numpy not %s'%str(type(coord)))

        if isinstance(coord, (float, int, complex)):
            coord = np.array(coord)

        if dim not in self.dims:
            self._dims.append(dim)
            self._coords.append(coord)
        else:
            raise ValueError('dim already in dims, cannot append to coords')

if __name__ == '__main__':

    coord = nddata_coord('x',slice(0,10,1))

    a = nddata_coord('a',slice(0,1,50e-3))
    b = np.r_[1:2:0.25]

#    d = nddata_coord_collection(a,coord)
    d = nddata_coord_collection(['a', 'x', 'b'],[a, coord, b])
    d.reorder(['x','a'])

    d.rename('a','test')

