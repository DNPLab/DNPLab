import numpy as _np

from .. import odnpData as _odnpData

from struct import unpack

headerSize = 32
blockHeaderSize = 28

header_fmt = '>llllllhhl'
blockHeader_fmt = '>hhhhlffff'


def importfid(path,filename):
    with open(path + filename,'rb') as f:
        headerString = f.read(headerSize)
        header = unpack(header_fmt,headerString)


        nblocks = header[0] # number of blocks in file
        ntraces = header[1] # number of traces per block
        npts = header[2] # number of elements per trace
        ebytes = header[3] # number of bytes per element
        tbytes = header[4] # number of bytes per traces
        bbytes = header[5] # number of bytes per block
        vers_id = header[6] # software version, file id status bits
        status = header[7] # status of whole file
        nbheaders =  header[8] # number of block headers per block

        # check if int or float
        isFloat = False
        if status & 0x08:
            is_float = True

        dataList = []
        for ix in range(nblocks):
            blockHeaderString = f.read(blockHeaderSize)
            blockHeader = unpack(blockHeader_fmt,blockHeaderString)

            scale = blockHeader[0] # scaling factor
            block_status = blockHeader[1] # status of data in block
            index = blockHeader[2] # block index
            mode = blockHeader[3] # mode of data in block
            ctcount = blockHeader[4] # ct value for FID
            lpval = blockHeader[5] # f2 (2D-f1) left phase in phasefile
            rpval = blockHeader[6] # f2 (2D-f1) right phase in phasefile
            lvl = blockHeader[7] # level drift correction
            tlt = blockHeader[8] # tilt drift correction

            blockDataString = f.read(tbytes)

            if isFloat:
                blockData = _np.array(unpack('>if'%(npts),blockDataString),dtype = complex)
            else:
                blockData = _np.array(unpack('>%ii'%(npts),blockDataString))
            data = blockData[0::2] + 1j*blockData[1:2]
            dataList.append(data)
        dataArray = _np.array(dataList).T

    return dataArray


def importProcpar(path,filename):
    paramDict = {}
    with open(path + filename,'r') as f:
        while True:
            line = f.readline()
            if line == '':
                return paramDict
            else:
                splitLine = line.rstrip().split(' ')

                varName = splitLine[0]

                firstValueLine = f.readline()
                valueLine = firstValueLine.rstrip().split(' ',1)
                numValues = valueLine[0]

                if numValues == '1':
                    try:
                        value = float(valueLine[1])
                    except:
                        value = valueLine[1]
                else:
                    # determine if value line is number
                    try:
                        valueLineNumberCheck = firstValueLine.rstrip().split(' ')
                        value = []
                        for eachValue in valueLineNumberCheck[1:]:
                            value.append(float(eachValue))

                    except:
                        value = [valueLine[1]]
                        for ix in range(int(numValues)-1):
                            otherValueLine = f.readline()
                            valueLine = otherValueLine.rstrip()
                            value.append(valueLine)

                finalLine = f.readline()
                enumValuesLine = finalLine.rstrip().split(' ')
                numEnumValues = enumValuesLine[0]
                if int(numEnumValues) > 0:
                   enumValues = enumValuesLine[1:]

                paramDict[varName] = value

def importVarian(path,filename,paramFilename = 'procpar'):
    '''
    '''

    paramDict = importProcpar(path,paramFilename)

    nmrFreq = paramDict['sfrq'] * 1.e6
    sw = paramDict['sw']
    npts = int(paramDict['np']/2)

    arraydim = paramDict['arraydim']

    dwellTime = 1./sw

    t = _np.r_[0.:int(npts)] * dwellTime
    
    data = importfid(path,filename)

    importantParamsDict = {}
    importantParamsDict['nmrFreq'] = nmrFreq
    print(_np.shape(data))
    print(len(t))
    output = _odnpData(data,[t],['t'],importantParamsDict)
    return output



