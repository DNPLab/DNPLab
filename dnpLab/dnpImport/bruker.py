import numpy as _np

from .. import dnpData as _dnpData

import os as _os

_dspfvs_table_10 = {
        2   : 44.7500,
        3   : 33.5000,
        4   : 66.6250,
        6   : 59.0833,
        8   : 68.5625,
        12  : 60.3750,
        16  : 69.5313,
        24  : 61.0208,
        32  : 70.0156,
        48  : 61.3438,
        64  : 70.2578,
        96  : 61.5052,
        128 : 70.3789,
        192 : 61.5859,
        256 : 70.4395,
        384 : 61.6263,
        512 : 70.4697,
        1024: 70.4849,
        1536: 61.6566,
        2048: 70.4924,
        }

_dspfvs_table_11 = {
        2   : 46.0000,
        3   : 36.5000,
        4   : 48.0000,
        6   : 50.1667,
        8   : 53.2500,
        12  : 69.5000,
        16  : 72.2500,
        24  : 70.1667,
        32  : 72.7500,
        48  : 70.5000,
        64  : 73.0000,
        96  : 70.6667,
        128 : 72.5000,
        192 : 71.3333,
        256 : 72.2500,
        384 : 71.6667,
        512 : 72.1250,
        1024: 72.0625,
        1536: 71.9167,
        2048: 72.0313,
        }

_dspfvs_table_12 = {
        2   : 46.311,
        3   : 36.530,
        4   : 47.870,
        6   : 50.229,
        8   : 53.289,
        12  : 69.551,
        16  : 71.600,
        24  : 70.184,
        32  : 72.138,
        48  : 70.528,
        64  : 72.348,
        96  : 70.700,
        128 : 72.524,
        192 : 71.3333,
        256 : 72.2500,
        384 : 71.6667,
        512 : 72.1250,
        1024: 72.0625,
        1536: 71.9167,
        2048: 72.0313,
        }

_dspfvs_table_13 = {
        2   : 2.750,
        3   : 2.833,
        4   : 2.875,
        6   : 2.917,
        8   : 2.938,
        12  : 2.958,
        16  : 2.969,
        24  : 2.979,
        32  : 2.984,
        48  : 2.989,
        64  : 2.992,
        96  : 2.995,
        }

def findGroupDelay(decim,dspfvs):
    '''
    '''
    groupDelay = 0
    if decim == 1:
        groupDelay = 0
    else:
        if dspfvs == 10:
            groupDelay = _dspfvs_table_10[int(decim)]
        elif dspfvs == 11:
            groupDelay = _dspfvs_table_11[int(decim)]
        elif dspfvs == 12:
            groupDelay = _dspfvs_table_12[int(decim)]
        elif dspfvs == 13:
            groupDelay = _dspfvs_table_13[int(decim)]
        else:
            print('dspfvs not defined')

    return groupDelay

def loadTitle(path, expNum = 1, titlePath = 'pdata/1',titleFilename = 'title'):
    '''
    Import Bruker Experiment Title File
    '''

    pathFilename = _os.path.join(path,str(expNum),titlePath,titleFilename) 

    with open(pathFilename,'r') as f:
        rawTitle = f.read()
    title = rawTitle.rstrip()

    return title

def loadAcqu(path, expNum = 1, paramFilename = 'acqus'):
    '''
    JCAMPDX file
    return dictionary of parameters
    '''

    pathFilename = path + str(expNum) + '/' + paramFilename

    # Import parameters
    with open(pathFilename,'r') as f:
        rawParams = f.read()

    # Split parameters by line
    lines = rawParams.strip('\n').split('\n')
    paramsDict = {}

    # Parse Parameters
    for line in lines:
        if line[0:3] == '##$':
            lineSplit = line[3:].split('= ')
            try:
                paramsDict[lineSplit[0]] = float(lineSplit[1])
            except:
                paramsDict[lineSplit[0]] = lineSplit[1]

    return paramsDict

def loadProc(path, expNum = 1, procNum = 1, paramFilename = 'procs'):
    '''
    '''

    pathFilename = path + str(expNum)  + '/pdata/' + str(procNum) + '/' + paramFilename

    # Import parameters
    with open(pathFilename,'r') as f:
        rawParams = f.read()

    # Split parameters by line
    lines = rawParams.strip('\n').split('\n')
    paramsDict = {}

    # Parse Parameters
    for line in lines:
        if line[0:3] == '##$':
            lineSplit = line[3:].split('= ')
            try:
                paramsDict[lineSplit[0]] = float(lineSplit[1])
            except:
                paramsDict[lineSplit[0]] = lineSplit[1]

    return paramsDict



def dirDataType(path,expNum):
    '''
    '''
    fullPath = path + '/' + str(expNum)

    dirList = _os.listdir(fullPath)

    if 'fid' in dirList:
        return 'fid'
    elif 'ser' in dirList: 
        if 'vdlist' in dirList:
            return 'ser'
        else:
            return 'serPhaseCycle'
    else:
        return ''

def importBruker(path,expNum,paramFilename = 'acqus'):
    '''
    '''
    dirType = dirDataType(path,expNum)

    if expNum is not None:
        fullPath = path + '/' + str(expNum)
    else:
        fullPath = path

    if dirType == 'fid':
        data = brukerFid(path,expNum,paramFilename)
        return data
    elif dirType == 'ser':
        data = brukerSer(path,expNum,paramFilename)
        return data
    elif dirType == 'serPhaseCycle':
        data = brukerSerPhaseCycle(path,expNum,paramFilename)
        return data
    else:
        raise ValueError
        Print('Could Not Identify Data Type in File')


def brukerFid(path,expNum,paramFilename = 'acqus'):
    '''
    '''
    paramsDict = loadAcqu(path, expNum, paramFilename)

    sw_h = paramsDict['SW_h'] # Spectral Width in Hz

    rg = paramsDict['RG'] # reciever gain
    decim = paramsDict['DECIM'] # Decimation factor of the digital filter
    dspfvs = paramsDict['DSPFVS'] # Digital signal processor firmware version
    bytorda = paramsDict['BYTORDA'] # 1 for big endian, 0 for little endian
    td = int(paramsDict['TD'])# points in time axes


    if bytorda == 0:
        endian = '<'
    else:
        endian = '>'

    raw = _np.fromfile(path + str(expNum) + '/fid',dtype = endian + 'i4')
    data = raw[0::2] + 1j * raw[1::2] # convert to complex

    groupDelay = findGroupDelay(decim,dspfvs)
    groupDelay = int(_np.ceil(groupDelay))

    t = 1./sw_h * _np.arange(0,int(td/2)-groupDelay)

    data = data[groupDelay:int(td/2)]

    data = data / rg

    importantParamsDict = {}
    importantParamsDict['nmrFreq'] = paramsDict['SFO1'] * 1e6
    output = _dnpData(data,[t],['t'],importantParamsDict)

    return output

def brukervdList(path,expNum):
    '''
    '''
    fullPath = path + str(expNum) + '/vdlist'

    with open(fullPath,'r') as f:
        raw = f.read()

#    lines = raw.strip('\n').split('\n')
    lines = raw.rstrip().rsplit()

    unitDict = {
            'n' : 1.e-9,
            'u' : 1.e-6,
            'm' : 1.e-3,
            'k' : 1.e3,
            }
    vdList = []
    for line in lines:
        if line[-1] in unitDict:
            value = float(line[0:-1]) * unitDict[line[-1]]
            vdList.append(value)
        else:
            value = float(line)
            vdList.append(value)


    vdList = _np.array(vdList)
    return vdList

def brukerSer(path,expNum,paramFilename = 'acqus'):
    '''
    '''
    paramsDict = loadAcqu(path, expNum, paramFilename)

    sw_h = paramsDict['SW_h'] # Spectral Width in Hz

    rg = paramsDict['RG'] # reciever gain
    decim = paramsDict['DECIM'] # Decimation factor of the digital filter
    dspfvs = paramsDict['DSPFVS'] # Digital signal processor firmware version
    bytorda = paramsDict['BYTORDA'] # 1 for big endian, 0 for little endian
    td = int(paramsDict['TD'])# points in time axes

    if bytorda == 0:
        endian = '<'
    else:
        endian = '>'

    raw = _np.fromfile(path + str(expNum) + '/ser',dtype = endian + 'i4')
    data = raw[0::2] + 1j * raw[1::2] # convert to complex

    groupDelay = findGroupDelay(decim,dspfvs)
    groupDelay = int(_np.ceil(groupDelay))

    t = 1./sw_h * _np.arange(0,int(td/2)-groupDelay)

    vdList = brukervdList(path,expNum)

    data = data.reshape(len(vdList),-1).T

    data = data[groupDelay:int(td/2),:]

    data = data / rg

    importantParamsDict = {}
    importantParamsDict['nmrFreq'] = paramsDict['SFO1'] * 1e6
    output = _dnpData(data,[t,vdList],['t','t1'],importantParamsDict)

    return output

def brukerSerPhaseCycle(path,expNum,paramFilename = 'acqus'):
    '''
    '''
    paramsDict = loadAcqu(path, expNum, paramFilename)

    sw_h = paramsDict['SW_h'] # Spectral Width in Hz

    rg = paramsDict['RG'] # reciever gain
    decim = paramsDict['DECIM'] # Decimation factor of the digital filter
    dspfvs = paramsDict['DSPFVS'] # Digital signal processor firmware version
    bytorda = paramsDict['BYTORDA'] # 1 for big endian, 0 for little endian
    td = int(paramsDict['TD'])# points in time axes

    if bytorda == 0:
        endian = '<'
    else:
        endian = '>'

    raw = _np.fromfile(path + str(expNum) + '/ser',dtype = endian + 'i4')
    data = raw[0::2] + 1j * raw[1::2] # convert to complex

    groupDelay = findGroupDelay(decim,dspfvs)
    groupDelay = int(_np.ceil(groupDelay))

    t = 1./sw_h * _np.arange(0,int(td/2)-groupDelay)


    length1d = int((_np.ceil(td/256.)*256)/2)
#    print length1d
    data = data.reshape(-1,int(length1d)).T

    data = data[groupDelay:int(td/2),:]

    # Assume phase cycle is 0, 90, 180, 270
    data = data[:,0] + 1j*data[:,1] - data[:,2] - 1j*data[:,3]
    data = data / rg

    importantParamsDict = {}
    importantParamsDict['nmrFreq'] = paramsDict['SFO1'] * 1e6

    output = _dnpData(data,[t],['t'],importantParamsDict)
    return output


def importBrukerDir(path):
    '''
    '''

    dirFiles = [x for x in _os.listdir(path) if _os.path.isdir(_os.path.join(path,x))]
    print(dirFiles)

    dataDict = {}
    for expNum in dirFiles:
        try:
            tempData = importBruker(path,expNum)
            dataDict[expNum] = tempData
        except:
            pass
#            print('%s is not a valid data directory'%(expNum))
    return dataDict




if __name__ == "__main__":
    pass
