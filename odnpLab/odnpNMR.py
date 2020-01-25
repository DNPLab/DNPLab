from . import odnpData as _odnpData
import numpy as _np


_defaultFourierTransformParameters = {
        'axes': 't',
        'zeroFillFactor' : 2, # zero fill axes to 
        'shift' : True, # shift axes to center zero frequency
        'convert2ppm' : True, # convert axes to ppm
        }

_defaultAlignParameters = {
        'axes': 't',
        }

_defaultWindowParameters = {
        'linewidth' : 10,
        'axes' : 't',
        }

_defaultIntegrateParameters = {
        'axes' : 't',
        'integrateCenter' : 0,
        'integrateWidth' : 100,
        }

def returnData(allData):
    '''
    Determine what the data type is for processing
    '''

    isDict = False
    if isinstance(allData,_odnpData):
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

def fourierTransform(allData, procParameters):
    '''
    Perform Fourier Transform down given dimension
    assumes dt = t[1] - t[0]

    Example:
    data.ft('t')
    '''

    # Determine if data is dictionary or odnpData object
    data, isDict = returnData(allData)

    requiredList = ['axes','zeroFillFactor','shift','convert2ppm']

    procParameters = updateParameters(procParameters,requiredList,_defaultFourierTransformParameters)

    axesLabel = procParameters['axes']
    zeroFillFactor = procParameters['zeroFillFactor']
    shift = procParameters['shift']
    convert2ppm = procParameters['convert2ppm']

    index = data.axesLabels.index(axesLabel)
    dt = data.axes[index][1] - data.axes[index][0]
    n_pts = zeroFillFactor*len(data.axes[index])
    f = (1./(n_pts*dt))*_np.r_[0:n_pts]
    if shift == True:
        f -= (1./(2*dt))

    if convert2ppm:
        nmrFrequency = data.params['nmrFreq']
        f /= (nmrFrequency / 1.e6)

    data.data = _np.fft.fft(data.data,n=n_pts,axis=index)
    if shift:
        data.data = _np.fft.fftshift(data.data,axes=index)
    data.axes[index] = f

    procStepString = 'ft:'
    for requiredParameter in requiredList:
        procStepString += requiredParameter + ',' + str(procParameters[requiredParameter]) + '\n'
    data.procList.append(procStepString)

    if isDict:
        allData['proc'] = data
        return allData
    else:
        return data

def window(allData,procParameters):
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

    data, isDict = returnData(allData)

    requiredList = ['axes','linewidth']

    procParameters = updateParameters(procParameters,requiredList,_defaultWindowParameters)

    axesLabel = procParameters['axes']
    linewidth = procParameters['linewidth']

    index = data.axesLabels.index(axesLabel)

    reshape_size = [1 for k in data.axesLabels]
    reshape_size[index] = len(data.axes[index])

#    if window_type == 'exp':
    window_array = _np.exp(-1.*data.axes[index]*linewidth).reshape(reshape_size)
        #NOTE Verify that this produces the correct linewidth, likely missing a factor of 2
#    else:
#        print('\'%s\' window type is not defined')
#        return
    window_array = _np.ones_like(data.data) * window_array
    data.data *= window_array

    procStepString = 'window:'
    for requiredParameter in requiredList:
        procStepString += requiredParameter + ',' + str(procParameters[requiredParameter]) + '\n'
    data.procList.append(procStepString)


    if isDict:
        allData['proc'] = data
        return allData
    else:
        return data

def integrate(allData,procParameters):
    '''
    Integrate data down given dimension
    
    Parameters:

    axes_label: str
        dimension to integrate
    integrateCenter:
        center for width of integration
    integrateWidth: 
        width of integration window

    Exampe:
    data.integrate(dataDict,procParameters)
    '''

    data, isDict = returnData(allData)

    requiredList = ['axes','integrateCenter','integrateWidth']

    procParameters = updateParameters(procParameters,requiredList,_defaultIntegrateParameters)

    axes = procParameters['axes']
    integrateCenter = procParameters['integrateCenter']
    integrateWidth = procParameters['integrateWidth']

    integrateMin = integrateCenter - _np.abs(integrateWidth)/2.
    integrateMax = integrateCenter + _np.abs(integrateWidth)/2.

    data = data.range(axes,integrateMin,integrateMax)

    data.sum(axes)

    if isDict:
        allData['proc'] = data
        return allData
    else:
        return data


def align(allData,procParameters):
    '''
    Alignment of NMR spectra down given axes dimension
    '''

    data, isDict = returnData(allData)

    if len(_np.shape(data.data)) != 2:
        print('Only 2-dimensional data supported')
        return

    requiredList = ['axes']

    procParameters = updateParameters(procParameters,requiredList,_defaultAlignParameters)

    alignAxesLabel = procParameters['axes']
    originalAxesOrder = data.axesLabels
    data.reorder(alignAxesLabel)
    axesIter = data.axesLabels[-1]

    refData = data[axesIter,0].data.reshape(-1)
    for ix in range(data.len(axesIter)):
        tempData = data[axesIter,ix].data.reshape(-1)

        corrData = _np.correlate(_np.abs(tempData),_np.abs(refData),mode='same')
        shiftIx = _np.argmax(corrData) - (len(corrData)/2) # subtract half length so spectrum is shifted relative to center, not edge
        shiftData = _np.roll(tempData,-1*shiftIx)
        data.data[:,ix] = shiftData
    data.reorder(originalAxesOrder)

    procStepString = 'align:'
    for requiredParameter in requiredList:
        procStepString += requiredParameter + ',' + str(procParameters[requiredParameter]) + '\n'
    data.procList.append(procStepString)

    if isDict:
        allData['proc'] = data
        return allData
    else:
        return data


def steps(allData):
    '''
    print processing steps
    '''

    if isinstance(allData,_odnpData):
        data = allData.copy()
    elif isinstance(allData,dict):
        isDict = True
        if 'proc' in allData:
            data = allData['proc'].copy()
        elif 'raw' in allData:
            data = allData['raw'].copy()
    else:
        print('Type not supported')
        return
    print('----------------')
    print('PROCESSING STEPS:')
    ix = 1
    for procStep in data.procList:
#        print('%i.)'%ix)
        procStep = procStep.split(':')
        print('----------------')
        print('%i.) '%ix + procStep[0] + ':')
        line = procStep[1].strip('\n').split('\n')
        for info in line:
            param_value = info.split(',')
            param = param_value[0]
            value = param_value[1]
            print(param + ', ' + value)
        ix += 1


