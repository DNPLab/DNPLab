'''dnpData object for storing N-dimensional data with coordinates
'''
import numpy as np
from collections.abc import MutableMapping
from copy import deepcopy

from .core import nddata

version = '1.0'

core_attrs_list = ['nmr_frequency']

class dnpData(nddata.nddata_core):
    '''dnpData Class for handling dnp data

    The dnpData class is inspired by pyspecdata nddata object which handles n-dimensional data, axes, and other relevant information together. 
    
    This class is designed to handle data and axes together so that performing NMR processing can be performed easily.

    Attributes:
    values (numpy.ndarray): Numpy Array containing data
    coords (list): List of numpy arrays containing axes of data
    dims (list): List of axes labels for data
    attrs (dict): Dictionary of parameters for data

    '''

    def __init__(self, values = np.r_[[]], coords = [], dims = [], attrs = {}, procList = []):
        '''dnpData Class __init__ method

        Args:
            data (numpy.ndarray): 
            coords (list): list of axes
            dims (list): list of strings which are names of axes
            attrs (dict): dictionary of parameters


        '''
        super().__init__(values, dims, coords, attrs)
        self.version = version
        self.proc_attrs = []
#
    def __repr__(self):
        '''Representation of dnpData object
        '''
        return 'nddata(values = {}, coords = {}, dims = {}, attrs = {})'.format(repr(self.values), repr(self.coords), repr(self.dims), repr(self.attrs))

    def __str__(self):
        '''String representation of dnpData object
        '''
        if len(self.attrs) < 20:
            return 'values:\n{}\ndims:\n{}\ncoords:\n{}\nattrs:\n{}\nproc_attrs:\n{}'.format(self.values, self.dims, self.coords, self.attrs,self.proc_attrs)
        else:
            core_attrs = {k:self.attrs[k] for k in core_attrs_list if k in core_attrs_list}
            num_additional_attrs = len(self.attrs) - len(core_attrs)
            return 'values:\n{}\ndims:\n{}\ncoords:\n{}\nattrs:\n{}\n + {} attrs'.format(self.values, self.dims, self.coords, core_attrs, num_additional_attrs)

    def add_proc_attrs(self,proc_attr_name,proc_dict):
        '''Stamp processing step to dnpdata object

        Args:
            proc_attr_name (str): Name of processing step (e.g. "fourier_transform"
            proc_dict (dict): Dictionary of processing parameters for this processing step.
        '''
        if not isinstance(proc_attr_name,str):
            raise ValueError('Processing step name must be string')
        if not isinstance(proc_dict,dict):
            raise ValueError('Processing dictionary must be dictionary')

        self.proc_attrs.append((proc_attr_name,proc_dict))

    def addAxes(self, dim, coord):
        '''Add new axesLabel to dnpData object with ax

        This function increases the dimension of the dnpData object by 1 with the axesValue parameter giving the axes

        Args:
            axesLabel (str): Name of new axis
            axesValue (float,int): Axes value for new dimension
        '''

        self.coords.append(dim, coord)
        self.values = np.expand_dims(self.values, -1)

    def autophase(self):
        '''Multiply dnpdata object by phase
        '''
        p = self.phase()
        self.values *= np.exp(-1j*p)
        if np.sum(np.real(self.values)) < 0:
            self.values *= -1.

    def concatenateAlong(self,newData,axesLabel):
        '''Concatenate new dnp data to original data along given axes label

        Args:
            newData (dnpData): data to be concatenated to dnp_data object
            axesLabel (str): axis to concatenate down
        '''
        reorderLabels = self.dims
        
        self.sort_dims()
        newData.sort_dims()


        if self.dims != newData.dims:
#            print 'ERROR' # NOTE determine how error handling will work
            raise ValueError('dims do not match')
#            return
        index = self.dims.index(axesLabel)

        self.values = np.concatenate((self.values,newData.values), axis = index)
        self.coords[axesLabel] = np.concatenate((np.array(self.coords[axesLabel]).reshape(-1),np.array(newData.coords[axesLabel]).reshape(-1)))

        self.reorder(reorderLabels)

    def phase(self,):
        '''Return phase of dnpData object

        Returns:
            phase (float,int): phase of data calculated from sum of imaginary divided by sum of real components
        '''
        return np.arctan(np.sum(np.imag(self.values))/np.sum(np.real(self.values)))

    def range(self,axesLabel,minValue,maxValue):
        '''Select range of data given axes values

        Args:
            axesLabel (str): Axes label for indexing
            minValue (float): Minimum axes value for indexing
            maxValue (float): Maximum axes value for indexing

        Returns:
            dnpData
        '''

        out = deepcopy(self)
        index = self.dims.index(axesLabel)

        min_list = list(out.coords[index] > minValue)
        max_list = list(out.coords[index] < maxValue)
        inRange = [min_list[i] and max_list[i] for i in range(len(min_list))]
        #NOTE if no values in range, will cause issues

        keep = [i for i, x in enumerate(inRange) if x]
        out.values = np.take(out.values,keep,axis=index)
        out.coords[index] = out.coords[index][keep]

        return out

    def squeeze(self):
        '''Remove all length 1 dimensions from data

        .. warning::
            Axes information is lost

        Example:
        data.squeeze()
        '''
        remove_axes = []
        for axes_ix,axes_value in enumerate(self.coords):
            if len(axes_value) == 1:
                remove_axes.append(axes_ix)
         
        reverse_remove_axes = remove_axes[::-1]
        for index_ix,index_value in enumerate(reverse_remove_axes):
            self.coords.pop(index_value)
            self.dims.pop(index_value)
            self.values = np.squeeze(self.values)

class dnpdata_collection(MutableMapping):
    '''Dictionary-like workspace object for storing dnpdata objects
    '''

    def __init__(self, *args, **kwargs):
        '''dnpdata_collection __init__ method

        Args:

        Example::

        '''
        self._processing_buffer = 'proc'

        self.__data_dict = {}

        if len(args) == 0:
            return
        elif len(args) == 1:
            if isinstance(args[0], dnpData):
                self.__data_dict['raw'] == dnpData
            elif isinstance(args[0], dict):
                data_dict = args[0]
                for key in data_dict:
                    if isinstance(data_dict[key], (dnpData, dict)):
                        self.__data_dict[key] = data_dict[key]
                    else:
                        raise TypeError('Each type in dict must be dnpData or dict')
            else:
                raise TypeError('Argument must be type dnpData')
        elif len(args) == 2:
            if isinstance(args[0], str) and isinstance(args[1], (dnpData, dict)):
                self.__data_dict[args[0]] = args[1]
            else:
                raise TypeError('If two arguments, first argument must be str and 2nd argument must be dnpData or dict')
        else:
            raise TypeError('Arguments not understood')

    def __getitem__(self, key):
        return self.__data_dict[key]

    def __setitem__(self, key, value):
        if (not isinstance(key, str)) or (not isinstance(value, (dict, dnpData))):
            raise TypeError('Key must be string and value must be dnpdata or dict')
        self.__data_dict[key] = value

    def __delitem__(self, key):
        del self.__data_dict[key]

    def __iter__(self):
        return iter(self.__data_dict)

    def __len__(self):
        return len(self.__data_dict)

    @property
    def processing_buffer(self):
        return self._processing_buffer

    @processing_buffer.setter
    def processing_buffer(self, new_processing_buffer):
        '''
        '''
        if isinstance(new_processing_buffer, str):
            self._processing_buffer = new_processing_buffer
        else:
            raise TypeError('Processing buffer must be type str, not %s'%str(type(new_processing_buffer)))

    def copy(self, a, b = None):
        '''Copy data
        '''

        if b is None:
            b = self.processing_buffer

        self[b] = self[a].copy()

    def move(self, a, b):
        '''Move data 
        '''

        self[b] = self.pop(a)

    def pop(self, b):
        return self.__data_dict.pop(b)

    def dict(self):
        return self.__data_dict

    def clear(self):
        '''
        '''
        self.__data_dict.clear()

    get = __getitem__

    def items(self):
        return self.__data_dict.items()

    def keys(self):
        return self.__data_dict.keys()

    def popitem(self):
        return self.__data_dict.popitem()

    def values(self):
        return self.__data_dict.values()

    def add(self, name, data):
        '''
        '''
        if (not isinstance(name, str)) or (not isinstance(data, (dnpData,dict))):
            raise TypeError('add takes two arguments, a string and dnpLab.odnpData type')
        self.__data_dict[name] = data

    def __repr__(self):
        return 'dnpdata_collection({})'.format(self.__data_dict)

    def __str__(self):
        return '{}\n'.format([(key,self[key].__str__()) for key in self.keys()])



def create_workspace(*args):
    return dnpdata_collection(*args)


if __name__ == '__main__':
    pass
