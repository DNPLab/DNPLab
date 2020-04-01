# Tim Keller
# Bridge12 Technologies, Inc
# Python Class for Handling ODNP Data
import numpy as np
#from matplotlib.pylab import *
from copy import deepcopy

version = '1.0'

#TODO:
# Fix Dictionary not being empty when odnp_data is initialized, problem with add_power function
# Force all axes to be nddata arrays, so that indexing is easier
# Add Alignment
# Add Fit down dimension

class odnpData:
    '''odnpData Class for handling odnp data

    The odnpData class is inspired by pyspecdata nddata object which handles n-dimensional data, axes, and other relevant information together. 
    
    This class is designed to handle data and axes together so that performing NMR processing can be performed easily.

    Attributes:
    data (numpy.ndarray): Numpy Array containing data
    axes (list): List of numpy arrays containing axes of data
    axesLabels (list): List of axes labels for data
    params (dict): Dictionary of parameters for data

    '''

    def __init__(self,data = np.r_[[]],axes = [],axesLabels = [],params = {},procList = []):
        '''odnpData Class __init__ method

        Args:
            data (numpy.ndarray): 
            axes (list): list of axes
            axesLabels (list): list of strings which are names of axes
            params (dict): dictionary of parameters


        '''

        self.data = data
        self.axes = axes
        self.axesLabels = axesLabels
        self.params = params

    def __add__(self,data):
        newData = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            newData.data += data
        elif isinstance(data,np.ndarray):
            newData.data = newData.data + data
        elif isinstance(data,odnpData):
            newData.data = newData.data + data.data
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return
        return newData

    def __radd__(self,data):
        return self.__add__(data)

    def __sub__(self,data):
        newData = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            newData.data -= data
        elif isinstance(data,np.ndarray):
            newData.data = newData.data - data
        elif isinstance(data,odnpData):
            newData.data = newData.data - data.data
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return
        return newData

    def __rsub__(self,data):
        newData = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            newData.data = data - newData.data
        elif isinstance(data,np.ndarray):
            newData.data = data - newData.data
        elif isinstance(data,odnpData):
            newData.data = data.data - newData.data
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return
        return newData


    def __mul__(self,data):
        newData = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            newData.data *= data
        elif isinstance(data,np.ndarray):
            newData.data = newData.data * data
        elif isinstance(data,odnpData):
            newData.data = newData.data * data.data
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return

        return newData

    def __rmul__(self,data):
        return self.__mul__(data)

    def __pow__(self,data):
        newData = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            newData.data = newData.data**data
        elif isinstance(data,np.ndarray):
            newData.data = newData.data**data
        elif isinstance(data,odnpData):
            newData.data = newData.data**data.data
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return
        return newData
        
    
    def __abs__(self):
        '''return absolute value of odnpData object
        Example::
           >> data = abs(data)
        '''
        out = deepcopy(self)
        out.data = np.abs(out.data)
        return out

    def real(self):
        out = deepcopy(self)
        out.data = np.real(out.data)
        return out

    def imag(self):
        out = deepcopy(self)
        out.data = np.imag(out.data)
        return out

    def abs(self):
        out = deepcopy(self)
        out.data = np.abs(out.data)
        return out
    def max(self):
        out = deepcopy(self)
        maxValue = np.max(out.data)
        return maxValue

    def __getitem__(self,index):
#        return self.data[index]
#        print index

        # Make sure format is correct
        if len(index) % 2 == 1:
            print('invalid indexing')
            return


        indexing_labels = index[0::2]
        indexing_list = index[1::2]

        indices = []

        for axes_ix, axes_label in enumerate(self.axesLabels):
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
                indices.append(slice(0,len(self.axes[axes_ix])))

        indices = tuple(indices)

        out = deepcopy(self)
        out.data = out.data[indices]
        for ix in range(len(out.axes)):
            out.axes[ix] = out.axes[ix][indices[ix]]


        return out

    def range(self,axesLabel,minValue,maxValue):

        out = deepcopy(self)
        index = self.axesLabels.index(axesLabel)

        min_list = list(out.axes[index] > minValue)
        max_list = list(out.axes[index] < maxValue)
        inRange = [min_list[i] and max_list[i] for i in range(len(min_list))]
        #NOTE if no values in range, will cause issues

        keep = [i for i, x in enumerate(inRange) if x]
        out.data = np.take(out.data,keep,axis=index)
        out.axes[index] = out.axes[index][keep]

        return out

    def __copy__(self):
        return deepcopy(self)
    def copy(self):
        return deepcopy(self)

    def __max__(self):
        '''
        '''
        #NOTE Not working
        out = deepcopy(self)
        out.data = np.max(out.data)
        return out
    def __len__(self):
        return np.size(self.data)

    def __str__(self):
        ''' String representation of odnpData object

        Returns:
            string (str): string representation of odnpData object
        '''
        string = 'Data Shape:' + repr(np.shape(self.data)) + '\n'
        string += 'Data Axes:' + repr(self.axesLabels) + '\n'
        string += 'Parameters:\n'# + repr(self.params)
        for key in self.params:
            if key != '*proc*':
                string += str(key) + ': ' + str(self.params[key]) + '\n'

        if '*proc*' in self.params:
            string += '*PROCESSING*\n'
            ix = 1
            for procStep in self.params['*proc*']:
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

    def __repr__(self):
        ''' Representation of odnpData object
        '''
        string = 'Data Shape:' + repr(np.shape(self.data)) + '\n'
        string += 'Data Axes:' + repr(self.axesLabels) + '\n'
        string += 'Parameters:\n'# + repr(self.params)
        for key in self.params:
            if key != '*proc*':
                string += str(key) + ': ' + str(self.params[key]) + '\n'

        if '*proc*' in self.params:
            string += '*PROCESSING*\n'
            ix = 1
            for procStep in self.params['*proc*']:
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

    def sort(self):
        '''Sort order of axes based on python list sorting for axes labels

        '''
        ix_sort = sorted(range(len(self.axesLabels)), key = lambda k: self.axesLabels[k])
        self.axes = [self.axes[ix] for ix in ix_sort]
        self.axesLabels = [self.axesLabels[ix] for ix in ix_sort]
        self.data = np.transpose(self.data,ix_sort)

    def reorder(self,axesLabels):
        '''
        Reorder array given a list of axes labels

        Args:
            axesLabels (list,tuple, str): Axes to reorder

        .. note::
            If not all axes are defined, they will be placed at the end of the axes labels in their original order
        '''
        if isinstance(axesLabels,str):
            axesLabels = [axesLabels]
        if isinstance(axesLabels,tuple):
            axesLabels = list(axesLabels)
        if not isinstance(axesLabels,list):
            print('axesLabels must be a list')
            return

        if sorted(axesLabels) != sorted(self.axesLabels):
            for label in axesLabels:
                if label not in self.axesLabels:
                    print('\'%s\' not in axes labels'%label)
                    return
            for label in self.axesLabels:
                if label not in axesLabels:
                    axesLabels.append(label)
        ix_reorder = [self.axesLabels.index(k) for k in axesLabels]
        self.axes = [self.axes[ix] for ix in ix_reorder]
        self.axesLabels = [self.axesLabels[ix] for ix in ix_reorder]
        self.data = np.transpose(self.data,ix_reorder)

    def rename(self, old_label, new_label):
        '''Rename axis
        
        Args:
            oldLabel (str): Axis label to be changed
            newLabel (str): New label for axes
        '''
        index = self.axesLabels.index(oldLabel)
        self.axesLabels[index] = newLabel

    def addAxes(self,axesLabel,axesValue):
        '''Add new axesLabel to odnpData object with ax

        This function increases the dimension of the odnpData object by 1 with the axesValue parameter giving the axes

        Args:
            axesLabel (str): Name of new axis
            axesValue (float,int): Axes value for new dimension
        '''
        if axesLabel in self.axesLabels:
            index = self.axesLabels.index(axesLabel)
            print('Axes %s already exists'%(str(axesLabel)))
        elif type(axesLabel) != str:
            index = axesLabel
            print('Axes label must be a string')
        else:
            self.axesLabels.append(axesLabel)
            self.axes.append(np.r_[axesValue])
            self.data = np.expand_dims(self.data,-1)

    def getAxes(self,axesLabel):
        '''Return given axes of odnpData object

        Args:
            axes_label (str): Axes to retrieve
        '''
        index = self.axesLabels.index(axesLabel)
        return self.axes[index]

    def index(self,axesLabel):
        '''Return index of given axes label

        Args:
            axesLabel (str): axis label to index
        '''
        return self.axesLabels.index(axesLabel)

    def concatenateAlong(self,newData,axesLabel):
        '''Concatenate new odnp data to original data along given axes label

        Args:
            newData (odnpData): data to be concatenated to odnp_data object
            axesLabel (str): axis to concatenate down
        '''
        reorderLabels = self.axesLabels
        
        self.sort()
        newData.sort()

        if self.axesLabels != newData.axesLabels:
#            print 'ERROR' # NOTE determine how error handling will work
            return
        index = self.axesLabels.index(axesLabel)

        self.data = np.concatenate((self.data,newData.data),axis = index)
        self.axes[index] = np.concatenate((self.axes[index],newData.axes[index]))

        self.reorder(reorderLabels)

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
#        index = self.axesLabels.index(axes_label)
#
#        plot_data = self.data
#
#        plot(self.axes[index],np.swapaxes(plot_data,0,index),*args,**kwargs)
#        xlabel(self.axesLabels[index])

    def squeeze(self):
        '''
        Remove all length 1 dimensions from data

        .. warning::
            Axes information is lost

        Example:
        data.squeeze()
        '''
        remove_axes = []
        for axes_ix,axes_value in enumerate(self.axes):
            if len(axes_value) == 1:
                remove_axes.append(axes_ix)
         
        reverse_remove_axes = remove_axes[::-1]
        for index_ix,index_value in enumerate(reverse_remove_axes):
            self.axes.pop(index_value)
            self.axesLabels.pop(index_value)
            self.data = np.squeeze(self.data)

    def sum(self,axesLabel):
        '''Perform sum down given dimension

        .. warning::
           Axis information is lost

        Args:
            axesLabel (str): Name of Axis to perform sum down

        .. code-block:: python

            data.sum('t')
        '''

        index = self.axesLabels.index(axesLabel)
        self.data = np.sum(self.data,axis = index)
        removedAxesLabel = self.axesLabels.pop(index)
        removedAxes = self.axes.pop(index)

    def autophase(self,):
        p = self.phase()
        self.data *= np.exp(-1j*p)
        if np.sum(np.real(self.data)) < 0:
            self.data *= -1.

    def phase(self,):
        '''Return phase of odnpData object

        Returns:
            phase (float,int): phase of data calculated from sum of imaginary divided by sum of real components
        '''
        return np.arctan(np.sum(np.imag(self.data))/np.sum(np.real(self.data)))

    def len(self,axesLabel):
        '''Return length of given dimension

        Args:
            axes_label (str): Axis to return length

        Example:
        data.len('t')
        '''
        index = self.axesLabels.index(axesLabel)

        return np.shape(self.data)[index]



if __name__ == '__main__':
    pass
