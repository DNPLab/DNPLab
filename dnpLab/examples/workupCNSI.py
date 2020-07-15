import numpy as np
from matplotlib.pylab import *
from scipy.io import loadmat, savemat
import odnpLab
from odnpLab.hydration import HydrationParameter, HydrationCalculator
from odnpLab.parameter import Parameter
import sys


# sys.path.append('/Users/thomascasey/odnplab')

class ProcParameter(Parameter):
    """Processing Parameters

    Attributes:
        eic (float): enhancement data integration window center
        eiw (float): enhancement data integration window width
        tic (float): T1 data integration window center
        tiw (float): T1 data integration window width
        verbose (bool): Whether verbose.

    """
    def __init__(self, init=None):  # TODO: enable manual adjustment
        super().__init__(init=init)
        eic = 0
        eiw = 100
        tic = 0
        tiw = 100
        self.eic, self.eiw, self.tic, self.tiw = eic, eiw, tic, tiw
        self.verbose = True


def process_cnsi(path: str, par: ProcParameter):
    """Process CNSI dataset

    Args:
        path (str): Path to CNSI data folder.
        par (ProcParameter): Processing parameters.

    Returns:
        Dictionary with T1, T1_power, E, E_power

    """

    eic, eiw, tic, tiw = par.eic, par.eiw, par.tic, par.tiw

    # Extract power settings from the experiment titles and later use them to make the power arrays
    Eplist = []
    for k in list(range(6, 27, 1)):
        title = odnpLab.odnpImport.bruker.loadTitle(path, expNum=k)
        splitTitle = title.split(' ')
        Eplist.append(float(splitTitle[-1]))

    T1plist = []
    for k in list(range(28, 33, 1)):
        title = odnpLab.odnpImport.bruker.loadTitle(path, expNum=k)
        splitTitle = title.split(' ')
        T1plist.append(float(splitTitle[-1]))

    addAttn = 38.5  # estimation of the attenuation between the directional coupler and power meter

    # convert .csv to .mat for older datasets
    try:
        power_raw = np.loadtxt(path + 'power.csv', delimiter=',', skiprows=1)

        time_list = power_raw[:, 0].reshape(-1, 1)
        power_list = power_raw[:, 1].reshape(-1, 1)

        power_dict = {'powerlist': power_list,
                      'timelist': time_list}

        # load t1_powers.csv
        t1_powers_raw = np.loadtxt(path + 't1_powers.csv', delimiter=',',
                                   skiprows=1)

        t1_time_list = t1_powers_raw[:, 0].reshape(-1, 1)
        t1_powers_list = t1_powers_raw[:, 1].reshape(-1, 1)

        t1_powers_dict = {'powerlist': t1_powers_list,
                          'timelist': t1_time_list}

        savemat(path + 'power.mat', power_dict)
        savemat(path + 't1_powers.mat', t1_powers_dict)

        print('Lists converted from csv to mat') if par.verbose else None
    except:
        print('Lists already in mat format') if par.verbose else None

    print('Started Ep processing...') if par.verbose else None
    EpNumList = list(range(5, 27, 1))
    print('ExpNumList = ', EpNumList) if par.verbose else None

    dataDir = odnpLab.odnpImport.bruker.importBrukerDir(path)

    power_t, power = odnpLab.odnpImport.power.importPower(path + 'power.mat')
    power_t, power = odnpLab.odnpImport.power.chopPower(power_t, power)
    expPowerList = power
    expPowerList = np.hstack(([0], expPowerList))

    data = odnpLab.odnpImport.power.assignPower(dataDir, EpNumList,
                                                expPowerList)

    dataDict = {'raw': data}

    dataDict = odnpLab.odnpNMR.removeOffset(dataDict, {})
    dataDict = odnpLab.odnpNMR.window(dataDict, {})
    dataDict = odnpLab.odnpNMR.fourierTransform(dataDict, {})

    # dataDict['proc'].data = dataDict['proc'].data * 1j

    phase = dataDict['proc']['t1', 1].phase()
    dataDict['proc'] *= np.exp(-1j * phase)

    dataDict = odnpLab.odnpNMR.integrate(dataDict, {'integrateCenter': eic,
                                                    'integrateWidth': eiw})

    # Normalize to first point
    dataDict['proc'].data /= dataDict['proc']['power', 0].data

    Ep_ = real(dataDict['proc'].data)
    Ep = Ep_[1:len(Ep_)]

    print('Ep processing Successful') if par.verbose else None

    Epows = np.multiply(-1, Eplist)
    Epows = np.add(Epows, addAttn)
    Epows = np.divide(Epows, 10)
    Epows = np.power(10, Epows)
    Epows = np.multiply((1e-3), Epows)

    print('Started T1 processing...') if par.verbose else None

    T1power_t, T1power = odnpLab.odnpImport.power.importPower(
        path + 't1_powers.mat')
    T1power_t, T1power = odnpLab.odnpImport.power.chopPower(T1power_t, T1power)

    T1pNumList = list(range(28, 33, 1))
    T1pNumList.append(304)
    print('T1Explist = ', T1pNumList) if par.verbose else None
    T1s = []
    T1pows = []
    for i in T1pNumList:

        data = odnpLab.odnpImport.bruker.importBruker(path, i)

        dataDict = {'raw': data}

        dataDict = odnpLab.odnpNMR.removeOffset(dataDict, {})
        dataDict = odnpLab.odnpNMR.window(dataDict, {})
        dataDict = odnpLab.odnpNMR.fourierTransform(dataDict, {})

        phase = dataDict['proc']['t1', 1].phase()
        dataDict['proc'] *= np.exp(-1j * phase)

        dataDict = odnpLab.odnpNMR.integrate(dataDict, {'integrateCenter': tic,
                                                        'integrateWidth': tiw})

        if dataDict['proc']['t1', -1].data < 0:
            dataDict['proc'] *= -1.

        dataDict = odnpLab.odnpFit.t1Fit(dataDict)

        T1s.append(dataDict['fit'].params['t1'])

    print('T1p processing Successful') if par.verbose else None

    T1pows = np.multiply(-1, T1plist)
    T1pows = np.add(T1pows, addAttn)
    T1pows = np.divide(T1pows, 10)
    T1pows = np.power(10, T1pows)
    T1pows = np.multiply((1e-3), T1pows)

    T1p = T1s[0:5]
    T10 = T1s[len(T1s) - 1]
    # T100 = 2.5

    return {'Epowers': np.array(Epows),
            'Ep': np.array(Ep),
            'T1powers': np.array(T1pows),
            'T1p': np.array(T1p),
            'T10': float(T10)}


if __name__ == '__main__':
    path = 'data/20190821_TW_4OH-TEMPO_500uM/'  # path to CNSI data folder
    ppar = ProcParameter()
    hpar = HydrationParameter()

    ppar.verbose = False

    rest = process_cnsi(path, ppar)
    t1, t1_power, e, e_power = rest['T1p'], rest['T1powers'], rest['Ep'], rest['Epowers']
    print(t1_power)

    hpar.slC = 500
    hpar.T10 = rest['T10']
    hpar.field = 380.4

    hc = HydrationCalculator(T1=t1,
                             T1_power=t1_power,
                             E=e,
                             E_power=e_power,
                             hp=hpar)
    print(hc.results)

# savemat('/Users/thomascasey/Documents/MATLAB/' + 'test.mat', {'odnp' : odnpData}, oned_as='column')
