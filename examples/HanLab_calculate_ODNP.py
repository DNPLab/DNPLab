"""ODNP Analysis for rb_dnp1 command at CNSI"""
### DO NOT EDIT
import datetime
import os
import re
import time

import numpy as np
from scipy.io import loadmat

import sys
sys.path.append('..')
import dnpLab as dnp
### DO NOT EDIT ABOVE ###

'''
INPUT YOUR PARAMS HERE:
'''
directory = '..data/topspin/' # path to data folder

Spin_Concentration = 100  # micro molar
Magnetic_Field = 348.5  # mT
T100 = 2.5  # T1 without spin label, without microwaves, commonly called T1,0(0)
smax_model = 'free' # ('tethered' or 'free')
t1_interp_method = 'linear' # ('linear' or 'second_order')


### DO NOT EDIT BELOW ###
def get_powers(directory,powerfile,ExpNums,bufferVal):
    """
    Args:
        path      (str)   :    A string of the odnp experiment folder. e.g. '../data/topspin/'
        powerfile (str)   :    A string specifying which file, either 'power' or 't1_powers'
        ExpNums   (list)  :    Folder numbers for experiments, either the enhancement folders list                     or the T1 folders list
        bufferVal (float) :    number that aligns something, 2.5 for 'power' and (20 * 2.5) for                        't1_power'
    Returns:
        power_list        :    array of power values taken from the .mat files
    """
    
    absTime = []
    for exp in ExpNums:
        opened = open(directory + str(exp) + os.sep + 'audita.txt')
        lines = opened.readlines()
        absStart = lines[8].split(' ')[2] + ' ' + lines[8].split(' ')[3]
        splitup = re.findall(r"[\w']+", absStart)
        absStart = datetime.datetime(int(splitup[0]), int(splitup[1]), int(splitup[2]), int(splitup[3]),
                                     int(splitup[4]), int(splitup[5]), int(splitup[6]))
        absStart = time.mktime(absStart.utctimetuple())  # this returns seconds since the epoch
        start = lines[8].split(' ')[3]
        start = start.split(':')  # hours,min,second
        hour = int(start[0], 10) * 3600
        minute = int(start[1], 10) * 60
        second = int(start[2].split('.')[0], 10)
        absStop = lines[6].split('<')[1].split('>')[0].split(' ')
        absStop = absStop[0] + ' ' + absStop[1]
        splitup = re.findall(r"[\w']+", absStop)
        absStop = datetime.datetime(int(splitup[0]), int(splitup[1]), int(splitup[2]), int(splitup[3]), int(splitup[4]),
                                    int(splitup[5]), int(splitup[6]))
        absStop = time.mktime(absStop.utctimetuple())  # this returns seconds since the epoch
        stop = lines[6].split(' ')[4]
        stop = stop.split(':')
        hour = int(stop[0], 10) * 3600
        minute = int(stop[1], 10) * 60
        second = int(stop[2].split('.')[0], 10)
        absTime.append((absStart, absStop))

    threshold = 20

    if os.path.isfile(directory + powerfile + '.mat'):  # This is a matlab file from cnsi
        print('Extracted powers from ' + powerfile + '.mat file')
        openfile = loadmat(directory + powerfile + '.mat')
        power = openfile.pop('powerlist')
        power = np.array([x for i in power for x in i])
        exptime = openfile.pop('timelist')
        exptime = np.array([x for i in exptime for x in i])
    elif os.path.isfile(directory + powerfile + '.csv'):  # This is a csv file
        print('Extracted powers from ' + powerfile + '.csv file')
        openfile = open(directory + powerfile + '.csv', 'r')
        lines = openfile.readlines()
        if len(lines) == 1:
            lines = lines[0].split('\r')  # this might not be what I want to do...
        lines.pop(0)
        timeList = []
        powerList = []
        for line in lines:
            exptime, power = line.split('\r')[0].split(',')
            timeList.append(float(exptime))
            powerList.append(float(power))
        exptime = np.array(timeList)
        power = np.array(powerList)

    #### Take the derivative of the power list
    step = exptime[1] - exptime[0]
    dp = []
    for i in range(len(power) - 1):
        dp.append((power[i + 1] - power[i]) / step)
    dp = abs(np.array(dp))
    ### Go through and threshold the powers
    timeBreak = []

    for i in range(len(dp)):
        if dp[i] >= threshold:
            timeBreak.append(exptime[i])
    timeBreak.sort()

    absTime.sort(key=lambda tup: tup[0])

    # align to the last spike
    offSet = absTime[-1][1] - timeBreak[-1] + bufferVal

    power_List = []
    for timeVals in absTime:
        start = int(timeVals[0] - offSet + bufferVal)
        stop = int(timeVals[1] - offSet - bufferVal)
        cutPower = []
        for k in range(0, len(exptime) - 1):
            if start <= exptime[k] <= stop:
                cutPower.append(power[k])
        powers = round(np.average(cutPower), 3)
        power_List.append(float(powers))
        
    return power_List


def find_peak(f, xmin, xmax):
    """ Find integer peak of f(x) for x \in [xmin, xmax]. T = O(log(n)).
    The frist derivative of f(x) must be monotonic. (by YL)

    +---------------+----------------+
    | xmax - xmin   | # calls of f   |
    +---------------+----------------+
    | 50            | 23             |
    +---------------+----------------+
    | 100           | 29             |
    +---------------+----------------+
    | 400           | 41             |
    +---------------+----------------+
    | 2000          | 53             |
    +---------------+----------------+

    Args:
        xmin: The minimum of possible x.
        xmax: The maximum of possible x
        f: Callable. Must return a number.
        *args: additional arguments for f.
        **kwargs: additional keyword arguments for f.

    Returns:
        A number x_ for which f(x_) >= f(x) for all x \in [xmin, xmax]

    """
    if xmin == xmax:
        return xmin, f(xmin)
    elif xmin > xmax:
        return xmax, f(xmax)

    ymin = f(xmin)
    ymax = f(xmax)

    if xmin == xmax - 1:
        return (xmin, ymin) if ymin >= ymax else (xmax, ymax)

    xmid = (xmin + xmax) // 2
    ymid = f(xmid)

    if ymin >= ymid >= ymax:
        xl, yl = find_peak(f, xmin+1, xmid-1)
        return (xl, yl) if yl >= ymin else (xmin, ymin)

    elif ymin <= ymid <= ymax:
        xr, yr = find_peak(f, xmid+1, xmax-1)
        return (xr, yr) if yr >= ymax else (xmax, ymax)

    else:
        x1, y1 = find_peak(f, xmin+1, xmid-1)
        x2, y2 = find_peak(f, xmid+1, xmax-1)
        (x1, y1) = (x1, y1) if y1 >= y2 else (x2, y2)
        return (x1, y1) if y1 >= ymid else (xmid, ymid)


def test_find_peak():
    # when xmin <= x_ <= xmax, return x_
    assert find_peak(lambda x:-x**2, -8, 8) == (0, 0)
    # when xmin <= xmax <= x_, return x_
    assert find_peak(lambda x:x, 0, 1) == (1, 1)
    # when f(a) == f(b) == f(x_) where a <= b, return a
    assert find_peak(lambda x:x**2, -1, 1) == (-1, 1)


def hanlab_calculate_odnp(directory:str, pars:dict, verbose=True):
    '''
    Args:
        directory: A string of the odnp experiment folder. e.g. '../data/topspin/'
        pars: A dictionary of the processing parameters (including integration width etc)
            Attr:
                ({integration_width  : int,   # set the default starting point to 20
                  spin_C             : float,
                  field              : float,
                  T100               : float,
                  smax_model         : str,   # ('tethered' or 'free')
                  t1_interp_method   : str    # ('linear' or 'second_order')
                })
        verbose: whether print intermediate outputs or not
    Returns:
        HydrationResults: A dnpLab.hydration.HydrationResults object.
    '''


    # folder number for p=0 point in ODNP set
    folder_p0 = 5

    # folder number for T1(0) in ODNP set
    folder_T10 = 304

    # list of folder numbers for Enhancement(p) points in ODNP set
    folders_Enhancements = range(6,27)

    # list of folder numbers for T1(p) points in ODNP set
    folders_T1s = range(28,33)

    total_folders = [folder_p0] + list(folders_Enhancements) + list(folders_T1s) + [folder_T10]

    T1 = []
    T1_stdd = []
    E = []
    for i, folder in enumerate(total_folders):

        data = dnp.dnpImport.topspin.import_topspin(directory, folder)
        workspace = dnp.create_workspace('raw',data)
        workspace.copy('raw','proc')

        dnp.dnpNMR.remove_offset(workspace,{})
        dnp.dnpNMR.window(workspace,{'linewidth' : 10})
        dnp.dnpNMR.fourier_transform(workspace,{'zero_fill_factor' : 2})

        if workspace['proc'].ndim == 2:
            workspace = dnp.dnpNMR.align(workspace, {})

        ## phase opt: optimize the phase
        curve = workspace['proc'].values

        phases = np.linspace(-np.pi / 2, np.pi / 2, 100).reshape(1, -1)
        rotated_data = (curve.reshape(-1, 1)) * np.exp(-1j * phases)
        success = (np.real(rotated_data) ** 2).sum(axis=0) / (
            (np.imag(rotated_data) ** 2).sum(axis=0))
        bestindex = np.argmax(success)

        phase = phases[0, bestindex]

        workspace['proc'] *= np.exp(-1j * phase)

        ## optCenter: find the optimized integration center
        def f_int(indx:int):
            y = sum(abs(dnp.dnpNMR.integrate(workspace['proc'], {'integrate_center' :  indx, 'integrate_width' : 10}).values))
            return y
        center, _ = find_peak(f_int, -50, 51)
        # intgrl_array = []
        # indxes = range(-50, 51)
        # workspace.copy('proc', 'proc0')
        # for indx in indxes:
        #     workspace = dnp.dnpNMR.integrate(workspace,{'integrate_center' :  indx, 'integrate_width' : 10})
        #     intgrl_array.append(sum(abs(workspace['proc'].values)))
        #     workspace.copy('proc0', 'proc')
        # cent = np.argmax(intgrl_array)
        # center = indxes[cent]

        workspace = dnp.dnpNMR.integrate(workspace,{'integrate_center' :  center, 'integrate_width' : pars['integration_width']})

        if  folder == folder_p0:
            p0 = np.real(workspace['proc'].values[0])
            print('Done with p0 folder...') if verbose else None
        elif  folder == folder_T10:
            workspace = dnp.dnpFit.t1Fit(workspace)
            T10 = workspace['fit'].attrs['t1']
            T10_stdd = workspace['fit'].attrs['t1_stdd']
            print('Done with T1(0) folder...') if verbose else None
        elif folder in folders_T1s:
            workspace = dnp.dnpFit.t1Fit(workspace)
            T1.append(workspace['fit'].attrs['t1'])
            T1_stdd.append(workspace['fit'].attrs['t1_stdd'])
            print('Done with T1(p) folder ' + str(folder) + '...') if verbose else None
        elif folder in folders_Enhancements:
            E.append(np.real(workspace['proc'].values[0]) / p0)
            print('Done with Enhancement(p) folder ' + str(folder) + '...') if verbose else None
        else:
            raise Exception('UnknownException')


    # get powers from .mat files
    # for E(p)
    powerfile = 'power'
    bufferVal = 2.5
    ExpNums = folders_Enhancements
    ## SEPARATE FUNCTION ##
    Epowers = get_powers(directory,powerfile,ExpNums,bufferVal)

    Epowers = np.add(Epowers, 21.9992)
    Epowers = np.divide(Epowers, 10)
    Epowers = np.power(10, Epowers)
    Epowers = np.multiply(1e-3, Epowers)

    E = np.array([Epowers, E])
    E = np.transpose(E)
    E = E[E[:,0].argsort()]

    Enhancement_powers = E[:,0]
    Enhancements = E[:,1]

    # for T1(p)
    powerfile = 't1_powers'
    bufferVal = 20 * 2.5
    ExpNums = folders_T1s
    ## SEPARATE FUNCTION ##
    T1powers = get_powers(directory,powerfile,ExpNums,bufferVal)

    T1powers = np.add(T1powers, 21.9992)
    T1powers = np.divide(T1powers, 10)
    T1powers = np.power(10, T1powers)
    T1powers = np.multiply(1e-3, T1powers)

    T1 = np.array([T1powers, T1, T1_stdd])
    T1 = np.transpose(T1)
    T1 = T1[T1[:,0].argsort()]

    T1_powers = T1[:,0]
    T1p = T1[:,1]
    T1p_stdd = T1[:,2]

    hydration = {
                 'E' : np.array(Enhancements),
                 'E_power' : np.array(Enhancement_powers),
                 'T1' : np.array(T1p),
                 'T1_power' : np.array(T1_powers)
                 }

    hydration.update({
                      'T10': T10,
                      'T100': pars['T100'],
                      'spin_C': pars['spin_C'],
                      'field': pars['field'],
                      'smax_model': pars['smax_model'],
                      't1_interp_method': pars['t1_interp_method']
                      })

    hydration_workspace = dnp.create_workspace()
    hydration_workspace.add('hydration_inputs', hydration)

    hydration_results = dnp.dnpHydration.hydration(hydration_workspace)

    hydration_results.update({'T1p_stdd': T1_stdd, 'T10_stdd': T10_stdd})

    return hydration_results


if __name__ == '__main__':
    pars = dict()
    pars.update({'integration_width'  : 20,
                 'spin_C'             : Spin_Concentration,
                 'field'              : Magnetic_Field,
                 'T100'               : T100,
                 'smax_model'         : smax_model,
                 't1_interp_method'   : t1_interp_method
               })

    hydration_results = hanlab_calculate_odnp(directory, pars)

    print('-----------------------')
    print('--------Results--------')
    print('krho = ' + str(round(hydration_results['krho'],2)) + ' (s-1M-1)')
    print('ksigma = ' + str(round(hydration_results['ksigma'],2)) + ' (s-1M-1), standard dev = ' + str(round(hydration_results['ksigma_stdd'],4)))
    print('coupling factor = ' + str(round(hydration_results['coupling_factor'],5)))
    print('tcorr = ' + str(round(hydration_results['tcorr'],2)) + ' ps')
    print('Dlocal = ' + str(round(hydration_results['Dlocal'] * 1e10,2)) + ' x 10^-10 (m^2/s)')
    print('-----------------------')
    print('-----------------------')

