import numpy as np

class nddata_axis(object):
    '''
    '''

    def __init__(self, dim, coord):
        '''

        Args:
            dim (str): name of dimension
        '''

        if isinstance(dim, str):
            self.__dim = dim
        else:
            raise ValueError('dim must be type str not %s'%(str(type(dim))))

        if isinstance(coord, np.ndarray):
            self.__coord = coord
        elif isinstance(coord, slice):
            self.__start = coord.start
            self.__stop = coord.stop
            self.__step = coord.step

        self.__size = self.array.size

    @property
    def start(self):
        '''
        '''
        return self.__start

    @property
    def stop(self):
        '''
        '''
        return self.__stop

    @property
    def step(self):
        '''
        '''
        return self.__step

    @property
    def array(self):
        '''Return axes as numpy array
        '''

        # Fast method
#        try:
#            return self.__array
#        except:
#            self.__array = np.r_[slice(self.start, self.stop, self.step)]
#            return self.__array
        # 2nd Attempt at fast method
        if hasattribute(self, '__array'):
            return self.__array
        except:
            self.__array = np.r_[slice(self.start, self.stop, self.step)]
            return self.__array

#        # Slow method, 56 times slower
#        return np.r_[slice(self.start, self.stop, self.step)]

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
        del self.__array


    @property
    def size(self):
        '''
        '''
        return self.array.size # 3 times faster when array is stored in object
#        return max(int((self.stop - self.start) / self.step),0) # Take max for case when direction is reversed


if __name__ == '__main__':

    coord = nddata_axis('x',slice(0,10,1))

    a = nddata_axis('x',slice(0,1,1e-3))

