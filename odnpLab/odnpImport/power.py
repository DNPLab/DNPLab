import numpy as _np

from scipy.io import loadmat as _loadmat

from .. import odnpData as _odnpData

def importPower(path,filename = ''):
    '''
    import powers file
    '''
    fullPath = path + filename

    if fullPath[-4:] == '.mat':
        rawDict = _loadmat(fullPath)
        t = rawDict['timelist'].reshape(-1)
        p = rawDict['powerlist'].reshape(-1)
    elif fullPath[-4:] == '.csv':
        raw = _np.loadtxt(fullPath,delimiter = ',',skiprows = 1)
        t = raw[:,0].reshape(-1)
        p = raw[:,1].reshape(-1)
    else:
        print('Could not identify power data type')
        return

    return t,p

def chopPower(t,p,threshold = 0.1):
    '''
    Use Derivative to chop Powers
    '''

    diffPower = _np.diff(p)

    step = [abs(x) > threshold for x in diffPower]

    correctedStep = []
    for ix in range(len(step) - 1):
        if step[ix] and step[ix+1]:
            correctedStep.append(False)
        elif step[ix] and not step[ix+1]:
            correctedStep.append(True)
        else:
            correctedStep.append(False)

    stepIndex = [0]
    for ix in range(len(correctedStep)):
        if correctedStep[ix]:
            stepIndex.append(ix)
    stepTupleList = []
    for ix in range(len(stepIndex)-1):
        stepTupleList.append((stepIndex[ix],stepIndex[ix+1]))

    averagePowerList = []
    averageTimeList = []
    for stepTuple in stepTupleList:
        averagePower = p[stepTuple[0]+1:stepTuple[1]]
        averagePower = _np.mean(averagePower)
        averagePowerList.append(averagePower)

        averageTime = (t[stepTuple[0]+1] + t[stepTuple[1]]) / 2.
        averageTimeList.append(averageTime)
    
    averagePowerArray = _np.array(averagePowerList)
    averageTimeArray = _np.array(averageTimeList)
    return averageTimeArray, averagePowerArray

def assignPower(dataDict,expNumList,powersList):
    '''
    Given a dictionary of odnpData objects with key being folder string,
    return the data with power values assigned to a new axis dimension
    '''


    doInitialize = True
    for ix,expNum in enumerate(expNumList):
        if str(expNum) in dataDict:
            if doInitialize:
                data = dataDict[str(expNum)]
                data.addAxes('power',powersList[ix])
                doInitialize = False
            else:
                tempData = dataDict[str(expNum)].copy()
                tempData.addAxes('power',powersList[ix])
                data.concatenateAlong(tempData,'power')


    return data

if __name__ == '__main__':
    from matplotlib.pylab import *
    #### Set Custom Matplotlib Parameters
    #matplotlib.rcParams['font.family'] = 'Myriad Pro' # font style, same as illustrator default
    matplotlib.rcParams['font.size'] = 24. # font size for axis

    matplotlib.rcParams['lines.linewidth'] = 1.
    matplotlib.rcParams['axes.linewidth'] = 2.
    matplotlib.rcParams['legend.fontsize'] = 14. # set legend font
    #matplotlib.rcParams['legend.fontsize'] = 20. # set legend font

    #matplotlib.rcParams['figure.subplot.bottom'] = 0.15
    matplotlib.rcParams['figure.subplot.bottom'] = 0.20
    matplotlib.rcParams['figure.subplot.top'] = .9
    #matplotlib.rcParams['figure.subplot.left'] = .125
    matplotlib.rcParams['figure.subplot.left'] = .20
    matplotlib.rcParams['figure.subplot.right'] = .9

    matplotlib.rcParams['xtick.major.size'] = 10 # sets tick thickness
    matplotlib.rcParams['xtick.major.width'] = 2 # sets tick thickness
    matplotlib.rcParams['xtick.minor.size'] = 5 # sets tick thickness
    matplotlib.rcParams['xtick.minor.width'] = 1 # sets tick thickness
    matplotlib.rcParams['xtick.direction'] = 'out' # sets tick thickness

    matplotlib.rcParams['ytick.major.size'] = 10 # sets tick thickness
    matplotlib.rcParams['ytick.major.width'] = 2 # sets tick thickness
    matplotlib.rcParams['ytick.minor.size'] = 5 # sets tick thickness
    matplotlib.rcParams['ytick.minor.width'] = 1 # sets tick thickness
    matplotlib.rcParams['ytick.direction'] = 'out' # sets tick thickness

    matplotlib.rcParams['pdf.fonttype'] = 42 # set font type so that I can edit with illustrator


    path = 'G:/My Drive/Exchange/Projects/0055 CPF-NIGM-0055 ODNP System/Software/Python/data/TEMPO_and_PEG_ODNP_data/TEMPO_and_PEG_ODNP_data/20191017_TW_4OHTEMPO_1p0mM/'
    filename = 'power.mat'

    t,p = importPower(path,filename)

    figure('raw powers')
    plot(t,p,linewidth = 2.)
    xlabel('Time (s)')
    averageTime,averagePower = chopPower(t,p)
#    figure('chopped Powers')
    plot(averageTime,averagePower,'bo',markersize = 6.)
    ylabel('Power (dB)')

    xlabel('Index')
    show()

