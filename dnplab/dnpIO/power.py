import numpy as _np
from scipy.io import loadmat as _loadmat

from .. import dnpData as _dnpData


def importPower(path, filename=""):
    """
    import powers file
    """
    fullPath = path + filename

    if fullPath[-4:] == ".mat":
        rawDict = _loadmat(fullPath)
        t = rawDict["timelist"].reshape(-1)
        p = rawDict["powerlist"].reshape(-1)

    elif fullPath[-4:] == ".csv":
        raw = _np.loadtxt(fullPath, delimiter=",", skiprows=1)
        t = raw[:, 0].reshape(-1)
        p = raw[:, 1].reshape(-1)

    else:
        print("Could not identify power data type")
        return

    return t, p


def chopPower(t, p, threshold=0.1):
    """
    Use Derivative to chop Powers
    """

    diffPower = _np.diff(p)

    step = [abs(x) > threshold for x in diffPower]

    correctedStep = []
    for ix in range(len(step) - 1):
        if step[ix] and step[ix + 1]:
            correctedStep.append(False)
        elif step[ix] and not step[ix + 1]:
            correctedStep.append(True)
        else:
            correctedStep.append(False)

    stepIndex = [0]
    for ix in range(len(correctedStep)):
        if correctedStep[ix]:
            stepIndex.append(ix)

    stepTupleList = []
    for ix in range(len(stepIndex) - 1):
        stepTupleList.append((stepIndex[ix], stepIndex[ix + 1]))

    averagePowerList = []
    averageTimeList = []
    for stepTuple in stepTupleList:
        averagePower = p[stepTuple[0] + 1 : stepTuple[1]]
        averagePower = _np.mean(averagePower)
        averagePowerList.append(averagePower)
        averageTime = (t[stepTuple[0] + 1] + t[stepTuple[1]]) / 2.0
        averageTimeList.append(averageTime)

    averagePowerArray = _np.array(averagePowerList)
    averageTimeArray = _np.array(averageTimeList)
    return averageTimeArray, averagePowerArray


def assignPower(dataDict, expNumList, powersList):
    """
    Given a dictionary of dnpData objects with key being folder string,
    return the data with power values assigned to a new axis dimension
    """

    doInitialize = True
    for ix, expNum in enumerate(expNumList):
        if str(expNum) in dataDict:
            if doInitialize:
                data = dataDict[str(expNum)]
                data.addAxes("power", powersList[ix])
                doInitialize = False
            else:
                tempData = dataDict[str(expNum)].copy()
                tempData.addAxes("power", powersList[ix])
                data.concatenateAlong(tempData, "power")

    return data
