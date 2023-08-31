import numpy as _np
import os

from .. import DNPData

from struct import unpack

headerSize = 32
blockHeaderSize = 28

header_fmt = ">llllllhhl"
blockHeader_fmt = ">hhhhlffff"


def array_coords(attrs):
    """Return array dimension coords from parameters dictionary

    Args:
        attrs (dict): Dictionary of procpar parameters

    Returns:
        tuple: dim and coord for array
    """

    try:
        array_delta = attrs["arraydelta"]
        array_dim = attrs["arraydim"]
        array_start = attrs["arraystart"]
        array_stop = attrs["arraystop"]
        array_value = attrs["array"]

        if array_dim != 1 and array_value != "":
            # OVJ has strange behaviour while saving... this can work but doesn't need to succeed thus the if clause
            try:
                coord = _np.array(attrs[array_value])
            except KeyError:
                if array_stop > array_start:
                    coord = _np.r_[array_start : array_stop + array_delta : array_delta]
                else:
                    array_max = attrs["arraymax"]
                    coord = _np.r_[array_start : array_max + array_delta : array_delta]
            dim = array_value
        elif array_dim != 1 and array_value == "":
            coord = _np.r_[array_start : array_stop + array_delta : array_delta]
            dim = "t1"
        else:
            coord = None
            dim = None

        return dim, coord

    except KeyError:
        coord = None
        dim = None

        return dim, coord


def import_fid(path, filename="fid"):
    """Import VnmrJ fid file

    Args:
        path (str): Directory of fid file
        filename (str): Name of fid file. "fid" by default

    Returns:
        numpy.ndarray: Array of data

    """
    with open(os.path.join(path, filename), "rb") as f:
        headerString = f.read(headerSize)
        header = unpack(header_fmt, headerString)

        nblocks = header[0]  # number of blocks in file
        ntraces = header[1]  # number of traces per block
        npts = header[2]  # number of elements per trace
        ebytes = header[3]  # number of bytes per element
        tbytes = header[4]  # number of bytes per traces
        bbytes = header[5]  # number of bytes per block
        vers_id = header[6]  # software version, file id status bits
        status = header[7]  # status of whole file
        nbheaders = header[8]  # number of block headers per block

        # check if int or float
        isFloat = False
        if status & 0x08:
            isFloat = True

        dataList = []
        for ix in range(nblocks):
            blockHeaderString = f.read(blockHeaderSize)
            blockHeader = unpack(blockHeader_fmt, blockHeaderString)

            scale = blockHeader[0]  # scaling factor
            block_status = blockHeader[1]  # status of data in block
            index = blockHeader[2]  # block index
            mode = blockHeader[3]  # mode of data in block
            ctcount = blockHeader[4]  # ct value for FID
            lpval = blockHeader[5]  # f2 (2D-f1) left phase in phasefile
            rpval = blockHeader[6]  # f2 (2D-f1) right phase in phasefile
            lvl = blockHeader[7]  # level drift correction
            tlt = blockHeader[8]  # tilt drift correction

            blockDataString = f.read(tbytes)

            if isFloat:
                blockData = _np.array(
                    unpack(">%if" % (npts), blockDataString), dtype=complex
                )

            else:
                blockData = _np.array(unpack(">%ii" % (npts), blockDataString))
            data = blockData[0::2] - 1j * blockData[1::2]  # minus sign for VNMRJ data
            dataList.append(data)
        dataArray = _np.array(dataList).T

    return dataArray


def import_procpar(path, filename="procpar"):
    """Import VnmrJ procpar parameters file

    Args:
        path (str): Directory of file

    Returns:
        dict: Dictionary of procpar parameters
    """
    paramDict = {}
    with open(os.path.join(path, filename), "r") as f:
        while True:
            line = f.readline()
            if line == "":
                return paramDict
            else:
                # Line 1: Name & Type Line
                splitLine = line.rstrip().split(" ")

                name = splitLine[0]
                subtype = splitLine[1]
                basictype = int(splitLine[2])
                maxvalue = float(splitLine[3])
                minvalue = float(splitLine[4])
                stepsize = float(splitLine[5])
                Ggroup = int(splitLine[6])
                Dgroup = int(splitLine[7])
                protection = int(splitLine[8])
                active = int(splitLine[9])
                intptr = int(splitLine[10])

                # 3 Cases:
                # basictype is 1 (real) -> Line 2 separated by spaces
                # basictype is 2 (string) & First number is 1 -> single string on same line inside double quotes
                # basic type is 2 (string) & first number is greater than 1 -> first element is on same line, subsequent elements are on next lines, strings are surrounded by double quotes

                # Line 2: Value line
                firstValueLine = f.readline()
                valueLine = firstValueLine.rstrip().split(" ")
                numValues = int(valueLine[0])

                if basictype == 1:
                    if numValues == 1:
                        value = float(valueLine[1])
                    else:
                        listFloats = []
                        for number in valueLine[1:]:
                            listFloats.append(float(number))

                        value = listFloats

                elif basictype == 2:
                    if numValues == 1:
                        value = valueLine[1].replace('"', "")
                    else:
                        listStrings = []
                        listStrings.append(valueLine[1].replace('"', ""))

                        for ix in range(numValues - 1):
                            nextValueLine = f.readline()
                            nextValue = nextValueLine.strip()
                            listStrings.append(nextValue.replace('"', ""))

                        value = listStrings

                finalLine = f.readline()
                enumValuesLine = finalLine.rstrip().split(" ")
                numEnumValues = enumValuesLine[0]

                if int(numEnumValues) == 1:
                    enumValues = enumValuesLine[1]
                else:
                    enumValues = enumValuesLine[1:]

                paramDict[name] = value


def import_vnmrj(path, fidFilename="fid", paramFilename="procpar"):
    """Import VnmrJ Data

    Args:
        path(str): path to experiment folder
        fidFilename(str): FID file name
        paramFilename(str): process parameter filename

    Returns:
        dnpdata: data in dnpdata object
    """

    attrs = import_procpar(path, paramFilename)

    nmr_frequency = attrs["H1reffrq"] * 1.0e6
    sw = attrs["sw"]
    npts = int(attrs["np"] / 2)

    attrs["nmr_frequency"] = nmr_frequency

    dim, coord = array_coords(attrs)

    dwellTime = 1.0 / sw

    t = _np.r_[0.0 : int(npts)] * dwellTime
    dims = ["t2"]
    coords = [t]

    data = import_fid(path, fidFilename)

    if coord is not None:
        dims.append(dim)
        coords.append(coord)
    else:
        data = data.reshape(-1)

    # Assign data/spectrum type
    attrs["experiment_type"] = "nmr_spectrum"

    return DNPData(data, dims, coords, attrs)
