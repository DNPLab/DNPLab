# Bridge12 Technologies, Inc
# Python Class for Handling ODNP Data
import numpy as np
#from matplotlib.pylab import *
from copy import deepcopy

version = '1.0'

#TODO:
# Fix Dictionary not being empty when dnp_data is initialized, problem with add_power function
# Force all axes to be nddata arrays, so that indexing is easier
# Add Alignment
# Add Fit down dimension

class dnpData:
    '''dnpData Class for handling dnp data

    The dnpData class is inspired by pyspecdata nddata object which handles n-dimensional data, axes, and other relevant information together. 
    
    This class is designed to handle data and axes together so that performing NMR processing can be performed easily.

    Attributes:
    values (numpy.ndarray): Numpy Array containing data
    coords (list): List of numpy arrays containing axes of data
    dims (list): List of axes labels for data
    attrs (dict): Dictionary of parameters for data

    '''

    def __abs__(self):
        '''return absolute value of dnpData object
        Example::
           >> data = abs(data)
        '''
        out = deepcopy(self)
        out.values = np.abs(out.values)
        return out

    def __add__(self,data):
        newData = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            newData.values += data
        elif isinstance(data,np.ndarray):
            newData.values = newData.values + data
        elif isinstance(data,dnpData):
            newData.values = newData.values + data.values
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return
        return newData

    def __copy__(self):
        return deepcopy(self)

    def __getitem__(self,index):
        '''Indexing dnpData function
        '''
#        return self.data[index]
#        print index

        # Make sure format is correct
        if len(index) % 2 == 1:
            print('invalid indexing')
            return

        indexing_labels = index[0::2]
        indexing_list = index[1::2]

        indices = []

        for axes_ix, axes_label in enumerate(self.dims):
            if axes_label in indexing_labels:
                ix = indexing_labels.index(axes_label)
#                indices.append(indexing_list[ix])
                if indexing_list[ix] == -1:
                    indices.append(slice(indexing_list[ix],None))
                elif isinstance(indexing_list[ix],int):
                    indices.append(slice(indexing_list[ix],indexing_list[ix]+1))
                else:
                    indices.append(indexing_list[ix])
            else:
                indices.append(slice(0,len(self.coords[axes_ix])))

        indices = tuple(indices)

        out = deepcopy(self)
        out.values = out.values[indices]
        for ix in range(len(out.coords)):
            out.coords[ix] = out.coords[ix][indices[ix]]

        return out

    def __init__(self,values = np.r_[[]],coords = [],dims = [],attrs = {},procList = []):
        '''dnpData Class __init__ method

        Args:
            data (numpy.ndarray): 
            coords (list): list of axes
            dims (list): list of strings which are names of axes
            attrs (dict): dictionary of parameters


        '''
        self.version = version

        self.values = values
        self.coords = coords
        self.dims = dims
        self.attrs = attrs

    def __len__(self):
        return np.size(self.values)

    def __mul__(self,data):
        newData = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            newData.values *= data
        elif isinstance(data,np.ndarray):
            newData.values = newData.values * data
        elif isinstance(data,dnpData):
            newData.values = newData.values * data.values
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return

        return newData

    def __max__(self):
        ''' Return maximum value of numpy array
        '''
        #NOTE Not working
        out = deepcopy(self)
        out.values = np.max(out.values)
        return out

    def __pow__(self,data):
        newData = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            newData.values = newData.values**data
        elif isinstance(data,np.ndarray):
            newData.values = newData.values**data
        elif isinstance(data,dnpData):
            newData.values = newData.values**data.values
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return
        return newData

    def __radd__(self,data):
        return self.__add__(data)

    def __repr__(self):
        ''' Representation of dnpData object
        '''
        string = 'Data Shape:' + repr(np.shape(self.values)) + '\n'
        string += 'Data Axes:' + repr(self.dims) + '\n'
        string += 'Parameters:\n'# + repr(self.attrs)
        for key in self.attrs:
            if key != '*proc*':
                string += str(key) + ': ' + str(self.attrs[key]) + '\n'

        if '*proc*' in self.attrs:
            string += '*PROCESSING*\n'
            ix = 1
            for procStep in self.attrs['*proc*']:
                procStep = procStep.split(':')
                string += '%i.) '%ix + procStep[0] + ':\n'
                line = procStep[1].strip('\n').split('\n')
                for info in line:
                    param_value = info.split(',')
                    param = param_value[0]
                    value = param_value[1]
                    string += param + ', ' + value + '\n'
                ix += 1
        return string

    def __rmul__(self,data):
        return self.__mul__(data)

    def __rsub__(self,data):
        newData = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            newData.values = data - newData.data
        elif isinstance(data,np.ndarray):
            newData.values = data - newData.values
        elif isinstance(data,dnpData):
            newData.values = data.values - newData.values
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return
        return newData

    def __str__(self):
        ''' String representation of dnpData object

        Returns:
            string (str): string representation of dnpData object
        '''
        string = 'Data Shape:' + repr(np.shape(self.values)) + '\n'
        string += 'Data Axes:' + repr(self.dims) + '\n'
        string += 'Parameters:\n'# + repr(self.attrs)
        for key in self.attrs:
            if key != '*proc*':
                string += str(key) + ': ' + str(self.attrs[key]) + '\n'

        if '*proc*' in self.attrs:
            string += '*PROCESSING*\n'
            ix = 1
            for procStep in self.attrs['*proc*']:
                procStep = procStep.split(':')
                string += '%i.) '%ix + procStep[0] + ':\n'
                line = procStep[1].strip('\n').split('\n')
                for info in line:
                    param_value = info.split(',')
                    param = param_value[0]
                    value = param_value[1]
                    string += param + ', ' + value + '\n'
                ix += 1
        return string

    def __sub__(self,data):
        newData = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            newData.values -= data
        elif isinstance(data,np.ndarray):
            newData.values = newData.values - data
        elif isinstance(data,dnpData):
            newData.values = newData.values - data.values
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return
        return newData
        
    def abs(self):
        '''Return absolute value of data
        '''
        out = deepcopy(self)
        out.values = np.abs(out.values)
        return out

    def addAxes(self,axesLabel,axesValue):
        '''Add new axesLabel to dnpData object with ax

        This function increases the dimension of the dnpData object by 1 with the axesValue parameter giving the axes

        Args:
            axesLabel (str): Name of new axis
            axesValue (float,int): Axes value for new dimension
        '''
        if axesLabel in self.dims:
            index = self.dims.index(axesLabel)
            print('Axes %s already exists'%(str(axesLabel)))
        elif type(axesLabel) != str:
            index = axesLabel
            print('Axes label must be a string')
        else:
            self.dims.append(axesLabel)
            self.coords.append(np.r_[axesValue])
            self.values = np.expand_dims(self.values,-1)

    def autophase(self,):
        '''Automatically phase data
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
        
        self.sort()
        newData.sort()

        if self.dims != newData.dims:
#            print 'ERROR' # NOTE determine how error handling will work
            return
        index = self.dims.index(axesLabel)

        self.values = np.concatenate((self.values,newData.values),axis = index)
        self.coords[index] = np.concatenate((self.coords[index],newData.coords[index]))

        self.reorder(reorderLabels)

    def copy(self):
        return deepcopy(self)

    def getAxes(self,axesLabel):
        '''Return given axes of dnpData object

        Args:
            axes_label (str): Axes to retrieve
        '''
        index = self.dims.index(axesLabel)
        return self.coords[index]

    def imag(self):
        '''Return imaginary part of data
        '''
        out = deepcopy(self)
        out.values = np.imag(out.values)
        return out

    def index(self,axesLabel):
        '''Return index of given axes label

        Args:
            axesLabel (str): axis label to index
        '''
        return self.dims.index(axesLabel)

    def len(self,axesLabel):
        '''Return length of given dimension

        Args:
            axes_label (str): Axis to return length

        Example:
        data.len('t')
        '''
        index = self.dims.index(axesLabel)

        return np.shape(self.values)[index]

    def max(self):
        '''Return maximum value of data
        '''
        out = deepcopy(self)
        maxValue = np.max(out.values)
        return maxValue

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

    def real(self):
        '''Return real part of data
        '''
        out = deepcopy(self)
        out.values = np.real(out.values)
        return out

    def reorder(self,dims):
        '''Reorder array given a list of axes labels

        Args:
            dims (list,tuple, str): Axes to reorder

        .. note::
            If not all axes are defined, they will be placed at the end of the axes labels in their original order
        '''
        if isinstance(dims,str):
            dims = [dims]
        if isinstance(dims,tuple):
            dims = list(dims)
        if not isinstance(dims,list):
            print('dims must be a list')
            return

        if sorted(dims) != sorted(self.dims):
            for label in dims:
                if label not in self.dims:
                    print('\'%s\' not in axes labels'%label)
                    return
            for label in self.dims:
                if label not in dims:
                    dims.append(label)
        ix_reorder = [self.dims.index(k) for k in dims]
        self.coords = [self.coords[ix] for ix in ix_reorder]
        self.dims = [self.dims[ix] for ix in ix_reorder]
        self.values = np.transpose(self.values,ix_reorder)

    def rename(self, oldLabel, newLabel):
        '''Rename axis
        
        Args:
            oldLabel (str): Axis label to be changed
            newLabel (str): New label for axes
        '''
        index = self.dims.index(oldLabel)
        self.dims[index] = newLabel

    def sort(self):
        '''Sort order of axes based on python list sorting for axes labels

        '''
        ix_sort = sorted(range(len(self.dims)), key = lambda k: self.dims[k])
        self.coords = [self.coords[ix] for ix in ix_sort]
        self.dims = [self.dims[ix] for ix in ix_sort]
        self.values = np.transpose(self.values,ix_sort)



#    def plot(self,axes_label,*args,**kwargs):
#        '''
#        plot data down given dimension
#
#        Parameters:
#        axes_label: str
#            axis to plot down (will be x-axis)
#        *args:
#            numpy args
#        **kwargs
#            numpy kwargs
#
#        NOTES:
#        Use show() to view figure
#
#        Example:
#        data.plot('t')
#        '''
#
#        index = self.dims.index(axes_label)
#
#        plot_data = self.values
#
#        plot(self.coords[index],np.swapaxes(plot_data,0,index),*args,**kwargs)
#        xlabel(self.dims[index])

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

    def sum(self,axesLabel):
        '''Perform sum down given dimension

        .. warning::
           Axis information is lost

        Args:
            axesLabel (str): Name of Axis to perform sum down

        .. code-block:: python

            data.sum('t')
        '''

        index = self.dims.index(axesLabel)
        self.values = np.sum(self.values,axis = index)
        removedAxesLabel = self.dims.pop(index)
        removedAxes = self.coords.pop(index)

if __name__ == '__main__':
    pass
