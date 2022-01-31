import numpy as np
import os
import re
import time
import datetime
from scipy.io import loadmat


def get_powers(path, power_file, experiment_list):
    """
    Split power readings files into array of power measurements equal in length to number of spectra in dataset

    Args:
        path (str) : Path to base folder containing power file
        power_file (str) : filename, "power" or "t1_powers"
        experiment_list (list) : list of folder numbers of experiments corresponding to power_file

    Returns:
        power_list (list) : list of power readings equal in length to experiment_list
    """

    if power_file == "power.mat" or power_file == "power.csv" or power_file == "power":
        power_file = "power"
    elif (
        power_file == "t1_powers.mat"
        or power_file == "t1_powers.csv"
        or power_file == "t1_powers"
    ):
        power_file = "t1_powers"
    else:
        raise TypeError("power file not recognized")

    if power_file == "power":
        buffer = 2.5
    elif power_file == "t1_powers":
        buffer = 50
    else:
        raise TypeError("invalid filename")

    try:
        expTime = []
        absTime = []
        for exp in experiment_list:
            opened = open(os.path.join(path, str(exp), "audita.txt"))
            lines = opened.readlines()
            absStart = lines[8].split(" ")[2] + " " + lines[8].split(" ")[3]
            splitup = re.findall(r"[\w']+", absStart)
            absStart = datetime.datetime(
                int(splitup[0]),
                int(splitup[1]),
                int(splitup[2]),
                int(splitup[3]),
                int(splitup[4]),
                int(splitup[5]),
                int(splitup[6]),
            )
            absStart = time.mktime(absStart.utctimetuple())
            start = lines[8].split(" ")[3]
            start = start.split(":")
            hour = int(start[0], 10) * 3600
            minute = int(start[1], 10) * 60
            second = int(start[2].split(".")[0], 10)
            start = second + minute + hour
            absStop = lines[6].split("<")[1].split(">")[0].split(" ")
            absStop = absStop[0] + " " + absStop[1]
            splitup = re.findall(r"[\w']+", absStop)
            absStop = datetime.datetime(
                int(splitup[0]),
                int(splitup[1]),
                int(splitup[2]),
                int(splitup[3]),
                int(splitup[4]),
                int(splitup[5]),
                int(splitup[6]),
            )
            absStop = time.mktime(absStop.utctimetuple())
            stop = lines[6].split(" ")[4]
            stop = stop.split(":")
            hour = int(stop[0], 10) * 3600
            minute = int(stop[1], 10) * 60
            second = int(stop[2].split(".")[0], 10)
            stop = second + minute + hour
            expTime.append(stop - start)
            absTime.append((absStart, absStop))

        threshold = 20

        if os.path.isfile(os.path.join(path, power_file + ".mat")):
            print("Extracted powers from " + power_file + ".mat file")
            openfile = loadmat(os.path.join(path, power_file + ".mat"))
            power = openfile.pop("powerlist")
            power = np.array([x for i in power for x in i])
            exptime = openfile.pop("timelist")
            exptime = np.array([x for i in exptime for x in i])
        elif os.path.isfile(os.path.join(path, power_file + ".csv")):
            print("Extracted powers from " + power_file + ".csv file")
            openfile = open(os.path.join(path, power_file + ".csv", "r"))
            lines = openfile.readlines()
            if len(lines) == 1:
                lines = lines[0].split("\r")
            lines.pop(0)
            timeList = []
            powerList = []
            for line in lines:
                exptime, power = line.split("\r")[0].split(",")
                timeList.append(float(exptime))
                powerList.append(float(power))
            exptime = np.array(timeList)
            power = np.array(powerList)

        step = exptime[1] - exptime[0]
        dp = []
        for i in range(len(power) - 1):
            dp.append((power[i + 1] - power[i]) / step)
        dp = abs(np.array(dp))

        timeBreak = []
        for i in range(len(dp)):
            if dp[i] >= threshold:
                timeBreak.append(exptime[i])

        timeBreak.sort()
        absTime.sort(key=lambda tup: tup[0])
        offSet = absTime[-1][1] - timeBreak[-1] + buffer

        power_list = []
        for timeVals in absTime:
            start = int(timeVals[0] - offSet + buffer)
            stop = int(timeVals[1] - offSet - buffer)
            cutPower = []
            for k in range(0, len(exptime) - 1):
                if start <= exptime[k] <= stop:
                    cutPower.append(power[k])
            powers = round(np.average(cutPower), 3)
            power_list.append(float(powers))
    except:
        raise ImportError("Unable to read the power file")

    return power_list
