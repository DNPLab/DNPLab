import numpy as _np
import os

from .. import dnpData as _dnpData

from struct import unpack

headerSize = 32
blockHeaderSize = 28

header_fmt = '>llllllhhl'
blockHeader_fmt = '>hhhhlffff'

def array_coords(attrs):
    '''Return array dimension coords from parameters dictionary

    Args:
        attrs (dict): Dictionary of procpar parameters

    Returns:
        tuple: dim and coord for array

    '''

    dim = attrs['array'].value

    if dim == '':
        dim = 'array'

    array_delta = attrs['arraydelta'].value


    array_max = attrs['arraymax'].value
    array_flip = attrs['arrayflip'].value

    array_start = attrs['arraystart'].value
    array_stop = attrs['arraystop'].value

    array_elements = attrs['arrayelemts'].value
    array_d_scale = attrs['arraydscale'].value
    array_dodc = attrs['arraydodc'].value

    coord = _np.r_[array_start:array_stop:array_delta]

    return dim, coord

def importfid(path,filename):
    with open(os.path.join(path, filename),'rb') as f:
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
            isFloat = True

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
                blockData = _np.array(unpack('>%if'%(npts),blockDataString),dtype = complex)

            else:
                blockData = _np.array(unpack('>%ii'%(npts),blockDataString))
            data = blockData[0::2] + 1j*blockData[1::2]
            dataList.append(data)
        dataArray = _np.array(dataList).T

    return dataArray


def importProcpar(path,filename):
    paramDict = {}
    with open(os.path.join(path, filename),'r') as f:
        while True:
            line = f.readline()
            if line == '':
                return paramDict
            else:
                # Line 1: Name & Type Line
                splitLine = line.rstrip().split(' ')

                varName = splitLine[0]
                variable = procparVariable(varName)

                variable.subtype = splitLine[1]
                variable.basictype = int(splitLine[2])
                variable.maxvalue = float(splitLine[3])
                variable.minvalue = float(splitLine[4])
                variable.stepsize = float(splitLine[5])
                variable.Ggroup = int(splitLine[6])
                variable.Dgroup = int(splitLine[7])
                variable.protection = int(splitLine[8])
                variable.active = int(splitLine[9])
                variable.intptr = int(splitLine[10])

                # 3 Cases:
                # basictype is 1 (real) -> Line 2 separated by spaces
                # basictype is 2 (string) & First number is 1 -> single string on same line inside double quotes
                # basic type is 2 (string) & first number is greater than 1 -> first element is on same line, subsequent elements are on next lines, strings are surrounded by double quotes

                # Line 2: Value line
                firstValueLine = f.readline()
                valueLine = firstValueLine.rstrip().split(' ')
                numValues = int(valueLine[0])

                variable.numValues = numValues

                if variable.basictype == 1:
                    if variable.numValues == 1:
                        variable.value = float(valueLine[1])
                    else:
                        listFloats = []
                        for number in valueLine[1:]:
                            listFloats.append(float(number))

                        variable.value = listFloats

                elif variable.basictype == 2:
                    if variable.numValues == 1:
                        variable.value = valueLine[1].replace('"','')
                    else:
                        listStrings = []
                        listStrings.append(valueLine[1].replace('"',''))

                        for ix in range(variable.numValues - 1):
                            nextValueLine = f.readline()
                            nextValue = nextValueLine.strip()
                            listStrings.append(nextValue.replace('"',''))

                        variable.value = listStrings





#                if numValues == '1':
#                    try:
#                        value = float(valueLine[1])
#                    except:
#                        value = valueLine[1]
#                else:
#                    # determine if value line is number
#                    try:
#                        valueLineNumberCheck = firstValueLine.rstrip().split(' ')
#                        value = []
#                        for eachValue in valueLineNumberCheck[1:]:
#                            value.append(float(eachValue))
#
#                    except:
#                        value = [valueLine[1]]
#                        for ix in range(int(numValues)-1):
#                            otherValueLine = f.readline()
#                            valueLine = otherValueLine.rstrip()
#                            value.append(valueLine)

                finalLine = f.readline()
                enumValuesLine = finalLine.rstrip().split(' ')
                numEnumValues = enumValuesLine[0]
                variable.numEnumValues = numEnumValues
                if int(numEnumValues) == 1:
                    variable.enumValues = enumValuesLine[1]
                else:
                    variable.enumValues = enumValuesLine[1:]

                paramDict[variable.name] = variable


def importVarian(path, fidFilename='fid', paramFilename ='procpar'):
    """

    Args:
        path(str): path to experiment folder
        fidFilename(str): FID file name
        paramFilename(str): process parameter filename

    Returns:
        dnpData: data

    """

    paramDict = importProcpar(path,paramFilename)

    nmr_frequency = paramDict['H1reffrq'].value*1.e6
    sw = paramDict['sw'].value
    npts = int(paramDict['np'].value/2)

#    arraydim = int(paramDict['arraydim'].value)
    dim, coord = array_coords(paramDict)

    dwellTime = 1./sw

    t = _np.r_[0.:int(npts)] * dwellTime
    
    data = importfid(path, fidFilename)

    if coord.size == 1:
        data = data.reshape(-1)
        output = _dnpData(data,[t],['t'],{})
    else:
#        data = data.T
#        data = data.reshape(-1,arraydim)
#        output = _dnpData(data,[t,_np.array(range(arraydim))],['t','x'],{})
        output = _dnpData(data,[t,coord],['t',dim],{})

    importantParamsDict = {}
    importantParamsDict['nmr_frequency'] = nmr_frequency
    output.attrs = importantParamsDict
#    print(_np.shape(data))
#    print(len(t))
    return output


class procparVariable():
    '''
    Storage container for procpar varibles
    Line 1:
    name, subtype, basictype, naxvalue, minvalue, stepsizei, Ggroup, Dgroup, protection, active, intptr
    '''
    def __init__(self,variableName):
        self.name = variableName

    def __repr__(self):
        string = str(self.name)
        string += ' ' + str(self.subtype)
        string += ' ' + str(self.basictype)
        string += ' ' + str(self.maxvalue)
        string += ' ' + str(self.minvalue)
        string += ' ' + str(self.stepsize)
        string += ' ' + str(self.Ggroup)
        string += ' ' + str(self.Dgroup)
        string += ' ' + str(self.protection)
        string += ' ' + str(self.active)
        string += ' ' + str(self.intptr) + '\n'

        string += str(self.numValues)
        if self.basictype == 1:
            if self.numValues == 1:
                string += ' ' + str(self.value) + '\n'
            else:
                for each in self.value:
                    string += ' ' + str(each)
                string += '\n'
        else:
            if self.numValues == 1:
                string += ' "' + self.value + '"\n'
            else:
                for value in self.value:
                    string += ' "' + value + '"\n'
        string += str(self.numEnumValues)

        if self.numEnumValues == 1:
            string += ' ' + str(self.enumValues)
        else:
            for enumValues in self.enumValues:
                string += ' ' + str(enumValues)

        return string
