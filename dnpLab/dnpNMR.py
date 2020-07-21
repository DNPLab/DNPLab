from . import dnpData as _dnpData
import numpy as _np


_defaultFourierTransformParameters = {
        'dim': 't',
        'zero_fill_factor' : 2, # zero fill dim to 
        'shift' : True, # shift dim to center zero frequency
        'convert_to_ppm' : True, # convert dim to ppm
        }

_defaultAlignParameters = {
        'dim': 't',
        }

_defaultWindowParameters = {
        'linewidth' : 10,
        'dim' : 't',
        }

_defaultIntegrateParameters = {
        'dim' : 't',
        'integrate_center' : 0,
        'integrate_width' : 100,
        }

_defaultRemoveOffsetParameters = {
        'dim' : 't',
        'offset_points' : 10,
        }

def returnData(allData):
    '''
    Determine what the data type is for processing
    '''

    isDict = False
    if isinstance(allData,_dnpData):
        data = allData.copy()
    elif isinstance(allData,dict):
        isDict = True
        if 'proc' in allData:
            data = allData['proc'].copy()
        elif 'raw' in allData:
            data = allData['raw'].copy()
    else:
        print('Type not supported')
        raise ValueError()
    dataType = 'dataDict'
    return data, isDict

def updateParameters(procParameters, requiredList, defaultParameters):
    '''
    For all parameters without default parameters, assign parameter value to default

    '''
    updatedProcParameters = procParameters
    for requiredParameter in requiredList:
        if not requiredParameter in updatedProcParameters:
            updatedProcParameters[requiredParameter] = defaultParameters[requiredParameter]
            print('Required parameter "%s" not given.\nSetting "%s" to default value of:'%(requiredParameter,requiredParameter))
            print(defaultParameters[requiredParameter])

    return updatedProcParameters

def procString(name,procParameters,requiredList):
    '''
    '''
    procStepString = name
    for requiredParameter in requiredList:
        procStepString += requiredParameter + ',' + str(procParameters[requiredParameter]) + '\n'

    return procStepString

def stampProcStep(data,procStepString):
    '''
    '''
    if '*proc*' in data.attrs:
        data.attrs['*proc*'].append(procStepString)
    else:
        data.attrs['*proc*'] = [procStepString]
    return data

def removeOffset(allData, procParameters):
    '''Remove DC offset from FID by averaging the last few data points and subtracting the average

    Args:
        allData (dnpData, dict): Data container for data
        procParameters (dict,procParam): Processing Parameters

    .. code-block:: python

       procParameters['dim'] = 't'
       procParameters['offset_points'] = 10

       outData = dnpLab.dnpNMR.removeOffset(allData,procParameters)

    Returns:
        allData (dnpData, dict)
    '''
    # Determine if data is dictionary or dnpData object
    data, isDict = returnData(allData)

    #requiredList = ['dim','offset_points']
    requiredList = _defaultRemoveOffsetParameters.keys()
    procParameters = updateParameters(procParameters,requiredList,_defaultRemoveOffsetParameters)
    dim = procParameters['dim']
    offset_points = int(procParameters['offset_points'])

    offsetData = data['t',-1*offset_points:].values
    offsetData = offsetData.reshape(-1)
    offset = _np.mean(offsetData)

    data -= offset

    procStepName = 'Remove Offset:'
    procStepString = procString(procStepName,procParameters,requiredList)
    data = stampProcStep(data,procStepString)

    if isDict:
        allData['proc'] = data
        return allData
    else:
        return data


def fourierTransform(allData, procParameters):
    '''Perform Fourier Transform down dim dimension given in procParameters

    .. Note::
        Assumes dt = t[1] - t[0]

    Args:
        allData (dnpData, dict): Data container
        procParameters (dict, procParam): Processing parameters

    Returns:
        allData (dnpData, dict): Processed data in container
        
    Example:

    .. code-block:: python
    
        procParameters['dim'] = 't'
        procParameters['zero_fill_factor'] = 2
        procParameters['shift'] = True
        procParameters['convert_to_ppm'] = True

        allData = dnpLab.dnpNMR.fourierTransform(allData, procParameters)
    '''

    # Determine if data is dictionary or dnpData object
    data, isDict = returnData(allData)

    requiredList = ['dim','zero_fill_factor','shift','convert_to_ppm']

    procParameters = updateParameters(procParameters,requiredList,_defaultFourierTransformParameters)

    dimLabel = procParameters['dim']
    zero_fill_factor = procParameters['zero_fill_factor']
    shift = procParameters['shift']
    convert_to_ppm = procParameters['convert_to_ppm']

    index = data.dims.index(dimLabel)
    dt = data.coords[index][1] - data.coords[index][0]
    n_pts = zero_fill_factor*len(data.coords[index])
    f = (1./(n_pts*dt))*_np.r_[0:n_pts]
    if shift == True:
        f -= (1./(2*dt))

    if convert_to_ppm:
        nmr_frequency = data.attrs['nmr_frequency']
        f /= (nmr_frequency / 1.e6)

    data.values = _np.fft.fft(data.values,n=n_pts,axis=index)
    if shift:
        data.values = _np.fft.fftshift(data.values,axes=index)
    data.coords[index] = f

    procStepName = 'Fourier Transform:'
    procStepString = procString(procStepName,procParameters,requiredList)
    data = stampProcStep(data,procStepString)

    if isDict:
        allData['proc'] = data
        return allData
    else:
        return data

def window(allData,procParameters):
    '''Apply Apodization to data down given dimension
    
    Args:
        allData (dnpData, dict): data container
        procParameters (dict, procParam): parameter values

    .. note::
        Axis units assumed to be seconds

    Exampe:

    .. code-block:: python

        procParameters = {
                'linewidth' : 10,
                'dim' : 't',
                }
        allData = dnpLab.dnpNMR.window(allData,procParameters)
        
    '''

    data, isDict = returnData(allData)

    requiredList = ['dim','linewidth']

    procParameters = updateParameters(procParameters,requiredList,_defaultWindowParameters)

    dimLabel = procParameters['dim']
    linewidth = procParameters['linewidth']

    index = data.dims.index(dimLabel)

    reshape_size = [1 for k in data.dims]
    reshape_size[index] = len(data.coords[index])

    # Must include factor of 2 in exponential to get correct linewidth ->
    window_array = _np.exp(-1.*data.coords[index]*2.*linewidth).reshape(reshape_size)
    window_array = _np.ones_like(data.values) * window_array
    data.values *= window_array

    procStepName = 'window:'
    procStepString = procString(procStepName,procParameters,requiredList)

    data = stampProcStep(data,procStepString)

    if isDict:
        allData['proc'] = data
        return allData
    else:
        return data

def integrate(allData,procParameters):
    '''Integrate data down given dimension
    
    Args:
        allData (dnpData,dict): Data container
        procParameters (dict, procParam): Processing Parameters
            dim_label: str
                dimension to integrate
            integrate_center:
                center for width of integration
            integrate_width: 
                width of integration window

    Returns:
        allData (dnpData,dict): Processed data

    Exampe:
    .. code-block:: python

        procParameters = {
                    'dim' : 't',
                    'integrate_center' : 0,
                    'integrate_width' : 100,
                    }
        dnpLab.dnpNMR.integrate(allData,procParameters)
    '''

    data, isDict = returnData(allData)

    requiredList = ['dim','integrate_center','integrate_width']

    procParameters = updateParameters(procParameters,requiredList,_defaultIntegrateParameters)

    dim = procParameters['dim']
    integrate_center = procParameters['integrate_center']
    integrate_width = procParameters['integrate_width']

    integrateMin = integrate_center - _np.abs(integrate_width)/2.
    integrateMax = integrate_center + _np.abs(integrate_width)/2.

    data = data.range(dim,integrateMin,integrateMax)

    data.sum(dim)

    procStepName = 'Integrate:'
    procStepString = procString(procStepName,procParameters,requiredList)
    data = stampProcStep(data,procStepString)

    if isDict:
        allData['proc'] = data
        return allData
    else:
        return data


def align(allData,procParameters):
    '''
    Alignment of NMR spectra down given dim dimension
    '''

    data, isDict = returnData(allData)

    if len(_np.shape(data.values)) != 2:
        print('Only 2-dimensional data supported')
        return

    requiredList = ['dim']

    procParameters = updateParameters(procParameters,requiredList,_defaultAlignParameters)

    alignAxesLabel = procParameters['dim']
    originalAxesOrder = data.dims
    data.reorder(alignAxesLabel)
    dimIter = data.dims[-1]

    refData = data[dimIter,0].data.reshape(-1)
    for ix in range(data.len(dimIter)):
        tempData = data[dimIter,ix].values.reshape(-1)

        corrData = _np.correlate(_np.abs(tempData),_np.abs(refData),mode='same')
        shiftIx = _np.argmax(corrData) - (len(corrData)/2) # subtract half length so spectrum is shifted relative to center, not edge
        shiftData = _np.roll(tempData,-1*int(shiftIx))
        data.values[:,ix] = shiftData
    data.reorder(originalAxesOrder)

    procStepName = 'Align:'
    procStepString = procString(procStepName,procParameters,requiredList)
    data = stampProcStep(data,procStepString)

    if isDict:
        allData['proc'] = data
        return allData
    else:
        return data


def steps(allData):
    '''
    print processing steps
    '''

    string = ''
    data,isDict = returnData(allData)
    string += '----------------\n'
    string += 'PROCESSING STEPS:\n'
    ix = 1
    if '*proc*' in data.attrs:
        for procStep in data.attrs['*proc*']:
#            string += '%i.)'%ix + '\n'
            procStep = procStep.split(':')
            string += '----------------\n'
            string += '%i.) '%ix + procStep[0] + ':\n'
            line = procStep[1].strip('\n').split('\n')
            for info in line:
                param_value = info.split(',')
                param = param_value[0]
                value = param_value[1]
                string += param + ', ' + value + '\n'
            ix += 1
        return string
    else:
        return 'No Processing Steps Found'


