from __future__ import division
import numpy as np
import operator
import functools


#allowed_domains = ['FT','IFT','LT','ILT','Wavelet','IWavelet']
allowed_domains = ['FT','IFT']


class nddata_axis(object):
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

        return nddata_axis(self.dim, slice(start, stop, step))
        
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
        return 'nddata_axis(\'{}\', {})'.format(self.dim,self.array)

    def __str__(self):
        return '\'{}\':{}'.format(self.dim,str(self.array))
    
    def __add__(self, b):
        start = self.start + b
        stop = self.stop + b
        step = self.step

        if hasattr(self, '_array'):
            del self.array

        return nddata_axis(self.dim, slice(start, stop, step))

    __radd__ = __add__

    def __sub__(self, b):
        start = self.start - b
        stop = self.stop - b
        step = self.step

        if hasattr(self, '_array'):
            del self.array

        return nddata_axis(self.dim, slice(start, stop, step))

    def __rsub__(self, b):
        # switch start and stop
        start = b - self.stop
        stop = b - self.start 
        step = self.step

        if hasattr(self, '_array'):
            del self.array

        return nddata_axis(self.dim, slice(start, stop, step))

    def __mul__(self, b):
        start = b * self.start
        stop = b * self.stop
        step = b * self.step

        if hasattr(self, '_array'):
            del self.array

        return nddata_axis(self.dim, slice(start, stop, step))

    __rmul__ = __mul__

    def __truediv__(self, b):
        start = self.start / b
        stop = self.stop / b
        step = self.step / b

        if hasattr(self, '_array'):
            del self.array

        return nddata_axis(self.dim, slice(start, stop, step))

    def __rtruediv__(self, b):
        start = b / self.start
        stop = b / self.stop
        step = b / self.step

        if hasattr(self, '_array'):
            del self.array

        return nddata_axis(self.dim, slice(start, stop, step))

    def __array__(self):
        return self.array

    def __len__(self):
        return len(self.array)

    def __matmul__(self, b):
        return self.start + self.step * b

#class nddata_axis_collection(MutableMapping):
class nddata_axis_collection(object):
    def __init__(self, *args, **kwargs):

        self.__dims = []
        self.__coords = []

        for arg in args:
            if isinstance(arg, nddata_axis):
                self.__dims.append(arg.dim)
                self.__coords.append(arg)
            elif isinstance(arg, np.ndarray):
                raise ValueError('No dim given for numpy ndarray')

        for kwarg in kwargs: #key, value pairs only in order for python >=3.6
            dim = kwarg
            coord = kwarg[dim]
            if isinstance(coord, nddata_axis):
                self.__dims.append(dim)
                arg.dim = dim
                self.__coords.append(coord)
            if isinstance(coord, np.ndarray):
                self.__dims.append(kwarg)
                self.__coords.append(nddata_axis(dim,coord))
            else:
                raise TypeError('kwarg not understood')

    def index(self, key):
        return self.__dims.index(key)

    def __getitem__(self, key):
        return self.coords[self.index(key)]

    def __setitem__(self, key, value):
        if not isinstance(key, str):
            raise TypeError('key must be type str not %s'%str(type(key)))
        if not isinstance(value, (nddata_axis, np.ndarray)):
            raise TypeError('argument must be type nddata_axis or numpy ndarray not %s'%str(type(value)))

        if isinstance(value, np.ndarray):
            value = nddata_axis(key, value)

        # if key already in dims, overwrite
        if key in self.dims:
            index = self.index(key)
            self.__coord[index] = value
        else:
            self.__dims.append(key)
            self.__coords.append(value)

    def __delitem__(self, key):
        index = self.index(key)
        del self.__dims[index]
        del self.__coords[index]

    @property
    def dims(self):
        return self.__dims
    @property
    def coords(self):
        return self.__coords

    def __repr__(self):
        return 'nddata_axis_collection({})'.format(self.coords)

    def __str__(self):
        return '{}'.format(self.coords)

    def __iter__(self):
        return iter(self.coords)

    def __len__(self):
        return len(self.coords)

    @property
    def shape(self):
        return tuple(len(k) for k in self.coords)

    @property
    def size(self):
        return functools.reduce(operator.mul, [len(k) for k in self.coords], 1)

    def reorder(self, dims):
        if len(dims) != len(self.dims):
            raise ValueError('Length of dims do not match')
        for dim in dims:
            if dim not in self.dims:
                raise ValueError('%s not in dims'%str(dim))

#        new_order = [dims.index(dim) for dim in self.dims if dims in self.dims]

#        new_order = [dims.index(dim) for dim in self.dims]
#        self.__dims = dims
#        self.__coords = [self.coords.index(dims[x]) if dims[x] in self.dims else np.r_[[]] for x in range(len(dims))]

#    def transpose(self, axis):

if __name__ == '__main__':

    coord = nddata_axis('x',slice(0,10,1))

    a = nddata_axis('a',slice(0,1,50e-3))

    d = nddata_axis_collection(a,coord)

