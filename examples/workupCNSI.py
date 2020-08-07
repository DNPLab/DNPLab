import numpy as np
from matplotlib.pylab import *
from scipy.io import loadmat, savemat
import dnpLab
from dnpLab import create_workspace
from dnpLab.dnpHydration import HydrationParameter, HydrationCalculator, Parameter
import sys


# sys.path.append('/Users/thomascasey/dnplab')

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
        title = dnpLab.dnpImport.topspin.load_title(path, expNum=k)
        splitTitle = title.split(' ')
        Eplist.append(float(splitTitle[-1]))

    T1plist = []
    for k in list(range(28, 33, 1)):
        title = dnpLab.dnpImport.topspin.load_title(path, expNum=k)
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

    dataDir = dnpLab.dnpImport.topspin.import_topspin_dir(path)

    power_t, power = dnpLab.dnpImport.power.importPower(path + 'power.mat')
    power_t, power = dnpLab.dnpImport.power.chopPower(power_t, power)
    expPowerList = power
    expPowerList = np.hstack(([0], expPowerList))

    data = dnpLab.dnpImport.power.assignPower(dataDir, EpNumList, expPowerList)

    dataDict = create_workspace()
    dataDict.add('raw', data)
    dataDict.copy('raw', 'proc')

    dataDict = dnpLab.dnpNMR.remove_offset(dataDict, {})
    dataDict = dnpLab.dnpNMR.window(dataDict, {})
    dataDict = dnpLab.dnpNMR.fourier_transform(dataDict, {})

    # dataDict['proc'].data = dataDict['proc'].data * 1j

    phase = dataDict['proc']['t2', 1].phase()
    dataDict['proc'] *= np.exp(-1j * phase)

    dataDict = dnpLab.dnpNMR.integrate(dataDict, {'integrate_center': eic,
                                                    'integrate_width': eiw})

    # Normalize to first point
    dataDict['proc'].values /= dataDict['proc']['power', 0].values

    Ep_ = real(dataDict['proc'].values)
    Ep = Ep_[1:len(Ep_)]

    print('Ep processing Successful') if par.verbose else None

    Epows = np.multiply(-1, Eplist)
    Epows = np.add(Epows, addAttn)
    Epows = np.divide(Epows, 10)
    Epows = np.power(10, Epows)
    Epows = np.multiply((1e-3), Epows)

    print('Started T1 processing...') if par.verbose else None

    # T1power_t, T1power = dnpLab.dnpImport.power.importPower(
    #     path + 't1_powers.mat')
    # T1power_t, T1power = dnpLab.dnpImport.power.chopPower(T1power_t, T1power)

    T1pNumList = list(range(28, 33, 1))
    T1pNumList.append(304)
    print('T1Explist = ', T1pNumList) if par.verbose else None
    T1s = []
    # T1pows = []
    for i in T1pNumList:

        data = dnpLab.dnpImport.topspin.import_topspin(path, i)

        dataDict = create_workspace()
        dataDict.add('raw', data)
        dataDict.copy('raw', 'proc')

        dataDict = dnpLab.dnpNMR.remove_offset(dataDict, {})
        dataDict = dnpLab.dnpNMR.window(dataDict, {})
        dataDict = dnpLab.dnpNMR.fourier_transform(dataDict, {})

        phase = dataDict['proc']['t1', 1].phase()
        dataDict['proc'] *= np.exp(-1j * phase)

        dataDict = dnpLab.dnpNMR.integrate(dataDict, {'integrate_center': tic,
                                                        'integrate_width': tiw})

        if dataDict['proc']['t1', -1].values < 0:
            dataDict['proc'] *= -1.

        dataDict = dnpLab.dnpFit.t1Fit(dataDict)

        T1s.append(dataDict['fit'].attrs['t1'])

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


def get_results(path:str, ppar:ProcParameter, hpar:HydrationParameter):
    """
    Returns:
        HydrationResults

    """
    rest = process_cnsi(path, ppar)
    t1, t1_power, e, e_power = rest['T1p'], rest['T1powers'], rest['Ep'], rest['Epowers']

    hpar.T10 = rest['T10']

    hc = HydrationCalculator(T1=t1,
                             T1_power=t1_power,
                             E=e,
                             E_power=e_power,
                             hp=hpar)
    hc.run()
    return hc.results


if __name__ == '__main__':
    path = 'data/TEMPO_500uM/'  # path to CNSI data folder
    ppar = ProcParameter()
    ppar.verbose = False
    ppar.tic, ppar.tiw, ppar.eic, ppar.eiw = 0, 250, 0, 250

    hpar = HydrationParameter()
    hpar.smax_model = 'free'
    hpar.t1_interp_method = 'linear'
    hpar.spin_C = 500
    hpar.field = 348.5



# savemat('/Users/thomascasey/Documents/MATLAB/' + 'test.mat', {'dnp' : dnpData}, oned_as='column')
