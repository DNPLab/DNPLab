# Tim Keller
# Bridge12 Technologies, Inc
# Python Class for Handling ODNP Data
import numpy as np
from matplotlib.pylab import *
from copy import deepcopy

version = '1.0'

#TODO:
# Fix Dictionary not being empty when odnp_data is initialized, problem with add_power function
# Force all axes to be nddata arrays, so that indexing is easier
# Add Alignment
# Add Fit down dimension

class odnpData:
    def __init__(self,data = np.r_[[]],axes = [],axesLabels = [],params = {},procList = []):
        self.data = data
        self.axes = axes
        self.axesLabels = axesLabels
        self.params = params

    def __add__(self,data):
        new_data = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            new_data.data += data
        elif isinstance(data,np.ndarray):
            new_data.data = new_data.data + data
        elif isinstance(data,odnpData):
            new_data.data = new_data.data + data.data
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return
        return new_data

    def __radd__(self,data):
        return self.__add__(data)

    def __sub__(self,data):
        new_data = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            new_data.data -= data
        elif isinstance(data,np.ndarray):
            new_data.data = new_data.data - data
        elif isinstance(data,odnpData):
            new_data.data = new_data.data - data.data
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return
        return new_data

    def __rsub__(self,data):
        new_data = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            new_data.data = data - new_data.data
        elif isinstance(data,np.ndarray):
            new_data.data = data - new_data.data
        elif isinstance(data,odnpData):
            new_data.data = data.data - new_data.data
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return
        return new_data


    def __mul__(self,data):
        new_data = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            new_data.data *= data
        elif isinstance(data,np.ndarray):
            new_data.data = new_data.data * data
        elif isinstance(data,odnpData):
            new_data.data = new_data.data * data.data
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return

        return new_data

    def __rmul__(self,data):
        return self.__mul__(data)

    def __pow__(self,data):
        new_data = deepcopy(self)
        if isinstance(data,(int,float,complex)) and not isinstance(data,bool):
            new_data.data = new_data.data**data
        elif isinstance(data,np.ndarray):
            new_data.data = new_data.data**data
        elif isinstance(data,odnpData):
            new_data.data = new_data.data**data.data
        else:
            print('Cannot add, type not supported:')
            print(type(data))
            return
        return new_data
        
    
    def __abs__(self):
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

    def range(self,axes_label,min_value,max_value):

        out = deepcopy(self)
        index = self.axesLabels.index(axes_label)

        min_list = list(out.axes[index] > min_value)
        max_list = list(out.axes[index] < max_value)
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
        '''
        Sort order of axes based on python list sorting for axes labels

        '''
        ix_sort = sorted(range(len(self.axesLabels)), key = lambda k: self.axesLabels[k])
        self.axes = [self.axes[ix] for ix in ix_sort]
        self.axesLabels = [self.axesLabels[ix] for ix in ix_sort]
        self.data = np.transpose(self.data,ix_sort)

    def reorder(self,axesLabels):
        '''
        Reorder array given a list of axes labels

        Parameters:
        axesLabels: list, tuple, str
            Axes to reorder

        NOTES:
        If not all axes are defined, they will be placed at the end of the axes labels in the original order
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
        '''
        Rename axis
        
        Parameters:
        old_label: str
            Axis label to be changed
        new_label: str
            New label for axes
        '''
        index = self.axesLabels.index(old_label)
        self.axesLabels[index] = new_label

    def add_axes(self,axes_label,axes_value):
        if axes_label in self.axesLabels:
            index = self.axesLabels.index(axes_label)
        elif type(axes_label) == int:
            index = axes_label
        self.axesLabels.append(axes_label)
        self.axes.append(np.r_[axes_value])
        self.data = self.data[:,np.newaxis]

    def get_axes(self,axes_label):
        '''
        Return given axes array

        Parameters:
        axes_label: str
            Axes to retrieve

        '''
        index = self.axesLabels.index(axes_label)
        return self.axes[index]

    def index(self,axes_label):
        '''
        return index of given axes label

        Parameters:
        axes_label: str
            axis label to index
        '''
        return self.axesLabels.index(axes_label)

    def concatenate(self,new_data):
        ''' 
        Add ODNP data to ODNP data
        Concatenate new odnp data to original data
        '''
        reorder_labels = self.axesLabels
        
        self.sort()
        new_data.sort()

        if self.axesLabels != new_data.axesLabels:
#            print 'ERROR' # NOTE determine how error handling will work
            return
        diff_axes = []
        for ix,axes in enumerate(self.axes):
            if not (np.array_equal(axes,new_data.axes[ix])):
                diff_axes.append(ix)

        if len(diff_axes) > 1:
            print('Only 1 dimension can be different')
            return
        if len(diff_axes) == 0:
#            print 'This function does not duplicate data'
            print('dimension to concatenate along is ambiguous, use concatenate_along')
            return

#        index = self.axesLabels.index(axes_label)
        index = diff_axes[0]

        self.data = np.concatenate((self.data,new_data.data),axis = index)
        self.axes[index] = np.concatenate((self.axes[index],new_data.axes[index]))

        self.reorder(reorder_labels)

    def concatenate_along(self,new_data,axes_label):
        ''' 
        Concatenate new odnp data to original data along given dimension

        Parameters:
        new_data: odnp_data object
            data to be concatenated to odnp_data object
        axes_label: str
            axis to concatenate down
        '''
        reorder_labels = self.axesLabels
        
        self.sort()
        new_data.sort()

        if self.axesLabels != new_data.axesLabels:
#            print 'ERROR' # NOTE determine how error handling will work
            return
        index = self.axesLabels.index(axes_label)

        self.data = np.concatenate((self.data,new_data.data),axis = index)
        self.axes[index] = np.concatenate((self.axes[index],new_data.axes[index]))

        self.reorder(reorder_labels)

    def add_power(self,new_data,power):
        '''
        Add Power to data ODNP data
        
        Parameters:
        new_data: odnp_data object
        power: float, int

        '''

        # if odnp_data is empty
        if (self.axes == []) and (self.axesLabels == []) and np.array_equal(self.data,np.r_[[]]):
            self.axes = new_data.axes
            self.axesLabels = new_data.axesLabels
            self.data = new_data.data
#            for key in new_data.params:
#                self.params[key] = new_data.params[key]
            self.params = new_data.params
            self.add_axes('power',power)
            return

        # This wouldn't work, it would only work for empty data
        if not 'power' in self.axesLabels:
            print('No Power Axes in Original Data')
            self.add_axes('power',power)

        if not 'power' in new_data.axesLabels:
            new_data.add_axes('power',power)
        self.concatenate_along(new_data,'power')

#    def add_t1(self,new_data,t1):
    def window(self,axes_label,linewidth = 10.,window_type = 'exp'):
        '''
        Apply Apodization to data down given dimension
        
        Parameters:

        axes_label: str
            dimension to apply apodization
        linewidth:
            Linewidth for exponential apodization, Hz
        window type: str
            Type of window function to apply

        NOTES:
        Axis units assumed to be seconds

        Exampe:
        data.window('t',linewidth = 20)
        '''
        index = self.axesLabels.index(axes_label)

        reshape_size = [1 for k in self.axesLabels]
        reshape_size[index] = len(self.axes[index])

        if window_type == 'exp':
            window_array = np.exp(-1.*self.axes[index]*linewidth).reshape(reshape_size)
        else:
            print('\'%s\' window type is not defined')
            return
        window_array = np.ones_like(data) * window_array
        self.data *= window_array

    def plot(self,axes_label,*args,**kwargs):
        '''
        plot data down given dimension

        Parameters:
        axes_label: str
            axis to plot down (will be x-axis)
        *args:
            numpy args
        **kwargs
            numpy kwargs

        NOTES:
        Use show() to view figure

        Example:
        data.plot('t')
        '''

        index = self.axesLabels.index(axes_label)

        plot_data = self.data

        plot(self.axes[index],np.swapaxes(plot_data,0,index),*args,**kwargs)
        xlabel(self.axesLabels[index])


    def zero_fill(self,axes_label,n_pts):
        '''
        Zero fill down given dimension
        Parameters:
        axes_label: str
            axis to zero fill down
        n_pts: int
            new length of zero filled axes

        NOTES:
        Assumes axes dimension is linear

        Example:
        data.zero_fill('t',16384)
        '''
        index = self.axesLabels.index(axes_label)

        dx = self.axes[index][1] - self.axes[index][0]
        x0 = self.axes[index][0]
        x = np.r_[0:dx*n_pts:1j*n_pts] + x0

        self.axes[index] = x
        pad = []
        for this_axes_label in self.axesLabels:
            if this_axes_label == axes_label:
                pad.append((0,n_pts - self.len(this_axes_label)))
            else:
                pad.append((0,0))
        pad = tuple(pad)
        self.data = np.pad(self.data,pad,mode = 'constant',constant_values = 0)


    def ft(self, axes_label,zero_fill_factor = 1, fftshift = True):
        '''
        Perform Fourier Transform down given dimension
        assumes dt = t[1] - t[0]

        Example:
        data.ft('t')
        '''

        index = self.axesLabels.index(axes_label)
        dt = self.axes[index][1] - self.axes[index][0]
        n_pts = zero_fill_factor*len(self.axes[index])
        f = (1./(n_pts*dt))*np.r_[0:n_pts]
        if fftshift == True:
            f -= (1./(2*dt))

        self.data = np.fft.fft(self.data,n=n_pts,axis=index)
        if fftshift:
            self.data = np.fft.fftshift(self.data,axes=index)
        self.axes[index] = f

    def ift(self, axes_label,zero_fill_factor = 1, fftshift = True):
        '''
        Perform Inverse Fourier Transform down given dimension
        assumes dt = t[1] - t[0]
        
        Example:
        data.ift('t')
        '''

        index = self.axesLabels.index(axes_label)
        dt = self.axes[index][1] - self.axes[index][0]
        n_pts = zero_fill_factor*len(self.axes[index])
        f = (1./(n_pts*dt))*np.r_[0:n_pts]
        if fftshift == True:
            f -= (1./(2*dt))

        self.data = np.fft.ifft(self.data,n=n_pts,axis=index)
        if fftshift:
            self.data = np.fft.fftshift(self.data,axes=index)
        self.axes[index] = f

    def squeeze(self):
        '''
        Remove all length 1 dimensions from data
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

    def sum(self,axes_label):
        '''
        Perform sum down given dimension

        Parameters:
        axes_label:
            axis label to sum down

        Example:
        data.sum('t')
        '''

        index = self.axesLabels.index(axes_label)
        self.data = np.sum(self.data,axis = index)
        removed_axes_label = self.axesLabels.pop(index)
        removed_axes = self.axes.pop(index)

    def autophase(self,):
        p = self.phase()
        self.data *= np.exp(-1j*p)
        if np.sum(np.real(self.data)) < 0:
            self.data *= -1.

    def phase(self,):
        return np.arctan(np.sum(np.imag(self.data))/np.sum(np.real(self.data)))

    def len(self,axes_label):
        '''
        Return length of given dimension

        Parameters:
        axes_label: str
            axis to return length

        Example:
        data.len('t')
        '''
        index = self.axesLabels.index(axes_label)

        return np.shape(self.data)[index]

    def align(self,axes_label,ref_data):
        '''
        Align spectra by maximizing cross correlation
        '''

        if len(ref_data.axesLabels) > 1:
            print('must be 1d reference')
            return

        axes_label = ref_data.axesLabels[0]
        reorder_array = self.axesLabels


        new_order = self.axesLabels
        new_order.insert(0,new_order.remove(axes_label))
        self.reorder(new_order)

        list(np.ones(len(self.axesLabels)))
#        new_shape = [1 for x in range(len(self.axesLabels))]
        new_shape[0] = len(self.axes)
        new_shape = [len(self.axes)]

        self.data.reshape()

        s = np.shape(data.data)

        toAlign = deepcopy(self)

        cross_corr = np.correlate(np.abs(spec),np.abs(spec_ref),mode='same')

        shift_ix = np.argmax(cross_corr) - (len(cross_corr)/2) # subtract half length so spectrum is shifted relative to center, not edge

        
        spec_shift = np.roll(spec,-1*shift_ix)

def test3d(std_noise = 0.):
    x = np.r_[0:100]
    y = np.r_[0:100]
    z = np.r_[0:100]

    noise = std_noise * np.random.randn(len(x),len(y),len(z))
    gauss = np.exp(-1.*(x-50)**2./(10.**2))
    gauss_3d = gauss.reshape(-1,1,1) * gauss.reshape(1,-1,1) * gauss.reshape(1,1,-1)
    gauss_3d += noise

    test_data = odnp_data(gauss_3d,[x,y,z],['x','y','z'])

    return test_data

        

def import_kea(path,filename = '',num = 1 ,verbose = False):
    params_dict = {}
    with open(path + filename + '/%i/'%num + 'acqu.par','r') as f:
        raw_params = f.read()

    raw_params = raw_params.strip().split('\n')
    for line in raw_params:
        if verbose:
            print(line)
        key_value = line.split(' = ')
        try:
            params_dict[key_value[0]] = float(key_value[1])
        except:
            params_dict[key_value[0]] = key_value[1]

    raw_data = np.loadtxt(path + filename + '/%i/'%num + 'data.csv', delimiter = ',')

    t = raw_data[:,0]
    t = t / 1.e6 # convert from us to s

    temp_data = raw_data[:,1] - 1j*raw_data[:,2]
    nmr_params = {}
    nmrFreq = float(params_dict['b1Freq'].strip('d'))*1e6 # convert to Hz

    nmr_params['nmrFreq'] = nmrFreq
    data = odnp_data(temp_data,[t],['t'],nmr_params)
    data.kea_params = params_dict
    return data


if __name__ == '__main__':
    # Experiment List
    exp_num_list = [1,24,33,37,39]

    # Convert to Array
    exp_num_array = np.array(exp_num_list)

    # Convert to power in dBm
    power_dBm_array = exp_num_array - 2.
    cable_loss = 1.5 # cable power loss in dB
    power_dBm_array -= cable_loss

    # Convert to Watts
    power_W_array = 10.**((power_dBm_array-30.)/10.)
    power_W_array[0] = 0 # Correct first value

    # Define Path and Import a few Test data sets
    path = './ref_data/10mM_TEMPOL_water_11-20-2019/'
    data = import_kea(path,num = 39)
    data.add_axes('power',0)
    data2 = import_kea(path,num = 1)
#    data2.add_axes('power',6e-3)

#    data.concatenate(data2)
#    data.concatenate_along(data2,'power')

    # Pre-allocate odnp_data class for enhancement profile
    data_power = odnpData()
    for ix in range(len(exp_num_list)):
        temp_data = import_kea(path,num = exp_num_list[ix])
        data_power.add_power(temp_data,power_W_array[ix])

    # Apply apodization
#    data_power.window('t',25) # Hz line broadening

    # Fourier Transform down time dimension
#    data_power.ft('t', zero_fill_factor = 2, fftshift = True)

   # Rename axis to frequency after Fourier Transform
#    data_power.rename_axes('t','f')

#    data_power.autophase()


#    figure()
#    data_power.plot('f')
#    show()


#    figure()
#    data_power['power',-1].plot('f')
#    show()
    # Integrate down frequency dimension
#    data_power.sum('f')


    # Plot data
#    figure()
#    data_power.plot('power','bo')

#    show()

