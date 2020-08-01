from . import dnpData, dnpdata_collection
import numpy as _np

_default_fourier_transform_parameters = {
        'dim': 't',
        'zero_fill_factor' : 2, # zero fill dim to 
        'shift' : True, # shift dim to center zero frequency
        'convert_to_ppm' : True, # convert dim to ppm
        }

_defaultAlign_parameters = {
        'dim': 't',
        }

_default_window_parameters = {
        'linewidth' : 10,
        'dim' : 't',
        }

_default_integrate_parameters = {
        'dim' : 't',
        'integrate_center' : 0,
        'integrate_width' : 100,
        }

_default_remove_offset_parameters = {
        'dim' : 't',
        'offset_points' : 10,
        }

_default_autophase_parameters = {
        'method': 'arctan',
        }

def return_data(all_data):
    '''
    Determine what the data type is for processing
    '''

    is_workspace = False
    if isinstance(all_data,dnpData):
        data = all_data.copy()
    elif isinstance(all_data, dict):
        raise ValueError('Type dict is not supported')
#        isDict = True
#        if 'proc' in all_data:
#            data = all_data['proc'].copy()
#        elif 'raw' in all_data:
#            data = all_data['raw'].copy()
    elif isinstance(all_data, dnpdata_collection):
        is_workspace = True
        if all_data.processing_buffer in all_data.keys():
            data = all_data[all_data.processing_buffer]
        else:
            raise ValueError('No data in processing buffer')
    else:
        raise ValueError('Data type not supported')

    return data, is_workspace

def update_parameters(proc_parameters, requiredList, default_parameters):
    '''
    For all parameters without default parameters, assign parameter value to default

    '''
    updatedProc_parameters = proc_parameters
    for requiredParameter in requiredList:
        if not requiredParameter in updatedProc_parameters:
            updatedProc_parameters[requiredParameter] = default_parameters[requiredParameter]
            print('Required parameter "%s" not given.\nSetting "%s" to default value of:'%(requiredParameter,requiredParameter))
            print(default_parameters[requiredParameter])

    return updatedProc_parameters

def remove_offset(all_data, proc_parameters):
    '''Remove DC offset from FID by averaging the last few data points and subtracting the average

    Args:
        all_data (dnpData, dict): Data container for data
        proc_parameters (dict,procParam): Processing _parameters

    .. code-block:: python

       proc_parameters['dim'] = 't'
       proc_parameters['offset_points'] = 10

       outData = dnpLab.dnpNMR.remove_offset(all_data,proc_parameters)

    Returns:
        all_data (dnpData, dict)
    '''
    # Determine if data is dictionary or dnpData object
    data, isDict = return_data(all_data)

    requiredList = _default_remove_offset_parameters.keys()
    proc_parameters = update_parameters(proc_parameters,requiredList,_default_remove_offset_parameters)
    dim = proc_parameters['dim']
    offset_points = int(proc_parameters['offset_points'])

    offsetData = data['t',-1*offset_points:].values
    offsetData = offsetData.reshape(-1)
    offset = _np.mean(offsetData)

    data -= offset

    proc_attr_name = 'remove_offset'
    proc_dict = {k:proc_parameters[k] for k in proc_parameters if k in requiredList}
    data.add_proc_attrs(proc_attr_name, proc_dict)

    if isDict:
        all_data[all_data.processing_buffer] = data
        return all_data
    else:
        return data


def fourier_transform(all_data, proc_parameters):
    '''Perform Fourier Transform down dim dimension given in proc_parameters

    .. Note::
        Assumes dt = t[1] - t[0]

    Args:
        all_data (dnpData, dict): Data container
        proc_parameters (dict, procParam): Processing parameters

    Returns:
        all_data (dnpData, dict): Processed data in container
        
    Example:

    .. code-block:: python
    
        proc_parameters['dim'] = 't'
        proc_parameters['zero_fill_factor'] = 2
        proc_parameters['shift'] = True
        proc_parameters['convert_to_ppm'] = True

        all_data = dnpLab.dnpNMR.fourier_transform(all_data, proc_parameters)
    '''

    # Determine if data is dictionary or dnpData object
    data, isDict = return_data(all_data)

    requiredList = _default_fourier_transform_parameters.keys()
    proc_parameters = update_parameters(proc_parameters,requiredList,_default_fourier_transform_parameters)

    dimLabel = proc_parameters['dim']
    zero_fill_factor = proc_parameters['zero_fill_factor']
    shift = proc_parameters['shift']
    convert_to_ppm = proc_parameters['convert_to_ppm']

    index = data.dims.index(dimLabel)
    dt = data.coords[dimLabel][1] - data.coords[dimLabel][0]
    n_pts = zero_fill_factor*len(data.coords[dimLabel])
    f = (1./(n_pts*dt))*_np.r_[0:n_pts]
    if shift == True:
        f -= (1./(2*dt))

    if convert_to_ppm:
        nmr_frequency = data.attrs['nmr_frequency']
        f /= (nmr_frequency / 1.e6)

    data.values = _np.fft.fft(data.values,n=n_pts,axis=index)
    if shift:
        data.values = _np.fft.fftshift(data.values,axes=index)
    data.coords[dimLabel] = f

    proc_attr_name = 'fourier_transform'
    proc_dict = {k:proc_parameters[k] for k in proc_parameters if k in requiredList}
    data.add_proc_attrs(proc_attr_name, proc_dict)

    if isDict:
        all_data[all_data.processing_buffer] = data
        return all_data
    else:
        return data

def window(all_data,proc_parameters):
    '''Apply Apodization to data down given dimension
    
    Args:
        all_data (dnpData, dict): data container
        proc_parameters (dict, procParam): parameter values

    .. note::
        Axis units assumed to be seconds

    Exampe:

    .. code-block:: python

        proc_parameters = {
                'linewidth' : 10,
                'dim' : 't',
                }
        all_data = dnpLab.dnpNMR.window(all_data,proc_parameters)
        
    '''

    data, isDict = return_data(all_data)

    requiredList = _default_window_parameters.keys()
    proc_parameters = update_parameters(proc_parameters,requiredList,_default_window_parameters)

    dimLabel = proc_parameters['dim']
    linewidth = proc_parameters['linewidth']

    index = data.dims.index(dimLabel)

    reshape_size = [1 for k in data.dims]
    reshape_size[index] = len(data.coords[dimLabel])

    # Must include factor of 2 in exponential to get correct linewidth ->
    window_array = _np.exp(-1.*data.coords[dimLabel]*2.*linewidth).reshape(reshape_size)
    window_array = _np.ones_like(data.values) * window_array
    data.values *= window_array

    proc_attr_name = 'window'
    proc_dict = {k:proc_parameters[k] for k in proc_parameters if k in requiredList}
    data.add_proc_attrs(proc_attr_name, proc_dict)

    if isDict:
        all_data[all_data.processing_buffer] = data
        return all_data
    else:
        return data

def integrate(all_data,proc_parameters):
    '''_integrate data down given dimension
    
    Args:
        all_data (dnpData,dict): Data container
        proc_parameters (dict, procParam): Processing Parameters
            dim_label: str
                dimension to integrate
            integrate_center:
                center for width of integration
            integrate_width: 
                width of integration window

    Returns:
        all_data (dnpData,dict): Processed data

    Exampe:
    .. code-block:: python

        proc_parameters = {
                    'dim' : 't',
                    'integrate_center' : 0,
                    'integrate_width' : 100,
                    }
        dnpLab.dnpNMR.integrate(all_data,proc_parameters)
    '''

    data, isDict = return_data(all_data)

    requiredList = _default_integrate_parameters.keys()
    proc_parameters = update_parameters(proc_parameters,requiredList,_default_integrate_parameters)

    dim = proc_parameters['dim']
    integrate_center = proc_parameters['integrate_center']
    integrate_width = proc_parameters['integrate_width']

    integrateMin = integrate_center - _np.abs(integrate_width)/2.
    integrateMax = integrate_center + _np.abs(integrate_width)/2.

#    data = data.range(dim,integrateMin,integrateMax)
    print('here')
    print(dim)
    data = data[dim,(integrateMin,integrateMax)]

    data = data.sum(dim)

    proc_attr_name = 'integrate'
    proc_dict = {k:proc_parameters[k] for k in proc_parameters if k in requiredList}
    data.add_proc_attrs(proc_attr_name, proc_dict)

    if isDict:
        all_data[all_data.processing_buffer] = data
        return all_data
    else:
        return data


def align(all_data,proc_parameters):
    '''
    Alignment of NMR spectra down given dim dimension
    '''

    data, isDict = return_data(all_data)

    if len(_np.shape(data.values)) != 2:
        print('Only 2-dimensional data supported')
        return

    requiredList = _defaultAlign_parameters.keys()
    proc_parameters = update_parameters(proc_parameters,requiredList,_defaultAlign_parameters)

    alignAxesLabel = proc_parameters['dim']
    originalAxesOrder = data.dims
    data.reorder([alignAxesLabel])
    dimIter = data.dims[-1]

    refData = data[dimIter,0].values.reshape(-1)
#    for ix in range(data.len(dimIter)):
    for ix in range(len(data.coords[dimIter])):
        tempData = data[dimIter,ix].values.reshape(-1)

        corrData = _np.correlate(_np.abs(tempData),_np.abs(refData),mode='same')
        shiftIx = _np.argmax(corrData) - (len(corrData)/2) # subtract half length so spectrum is shifted relative to center, not edge
        shiftData = _np.roll(tempData,-1*int(shiftIx))
        data.values[:,ix] = shiftData
    data.reorder(originalAxesOrder)

    proc_attr_name = 'align'
    proc_dict = {k:proc_parameters[k] for k in proc_parameters if k in requiredList}
    data.add_proc_attrs(proc_attr_name, proc_dict)

    if isDict:
        all_data[all_data.processing_buffer] = data
        return all_data
    else:
        return data

def autophase(workspace, parameters):
    '''
    '''

    requiredList = _default_autophase_parameters.keys()
    parameters = update_parameters(parameters, requiredList, _default_autophase_parameters)

    data, is_workspace = return_data(workspace)
    phase = _np.arctan(_np.sum(_np.imag(data.values))/_np.sum(_np.real(data.values)))

    data.values *= _np.exp(-1j*phase)
    if _np.sum(_np.real(data.values)) < 0:
        data.values *= -1.

    proc_attr_name = 'autophase'
    proc_attrs = {k:parameters[k] for k in parameters if k in requiredList}
    data.add_proc_attrs(proc_attr_name, proc_attrs)

    if is_workspace:
        workspace[workspace.processing_buffer] = data
        return workspace
    else:
        return data
