"""Example processing Bruker ODNP data and calcualting hydration parameters"""

import sys
import os
import copy
import numpy as np

sys.path.append('..')
import dnpLab as dnp

# directory of base folder containing the numbered folders listed below
directory = '..topspin'

# folder number for p=0 point in ODNP set
folder_p0 = 5

# folder number for T1(0) in ODNP set
folder_T10 = 304

# list of folder numbers for Enhancement(p) points in ODNP set
folders_Enhancements = range(6,27)

# list of folder numbers for T1(p) points in ODNP set
folders_T1s = range(28,33)

# give input parameters and set options for dnpHydration
spin_C = 100
field = 348.5
T100 = 2.5
smax_model = 'tethered'
t1_interp_method = 'second_order'

# array of power readings for Enhancement(p) points
Enhancement_powers = np.array([0.0006454923080882520,
0.004277023425898170,
0.004719543572446050,
0.00909714298712173,
0.01344187403986090,
0.01896059941058610,
0.02101937603827090,
0.022335737104727900,
0.026029715703921800,
0.02917012237740640,
0.0338523245243911,
0.03820738749745440,
0.04733370907740660,
0.05269608016472140,
0.053790874615060400,
0.05697639350179900,
0.06435487925718170,
0.07909179437004270,
0.08958910066880800,
0.1051813598911370,
0.11617812912435900])

# array of power readings for T1(p) points
T1_powers = np.array([0.000589495934876689,
0.024242327290569100,
0.054429505156431400,
0.0862844940360515,
0.11617812912435900])

proc_params = dict()
proc_params['zero_fill_factor'] = 2
proc_params['window_linewidth'] = 10

int_params = dict()
int_params['integrate_center'] = 0
int_params['integrate_width'] = 20

save_name = 'ODNPdata_dnpHydrationResults'


#### Do not change the code below ####
print('Working...')

def workupPhaseOpt(workspace):

    curve = workspace['proc'].values
    
    phases = np.linspace(-np.pi / 2, np.pi / 2, 100).reshape(1, -1)
    rotated_data = (curve.reshape(-1, 1)) * np.exp(-1j * phases)
    success = (np.real(rotated_data) ** 2).sum(axis=0) / (
        (np.imag(rotated_data) ** 2).sum(axis=0))
    bestindex = np.argmax(success)

    return phases[0, bestindex]

def optCenter(workspace):

    intgrl_array = []
    indx = range(-50, 51)
    int_params = {'integrate_width': 10}
    for k in range(0, len(indx)):
        iterativeopt_workspace = copy.deepcopy(workspace)
        iterativeopt_workspace = dnp.dnpNMR.integrate(iterativeopt_workspace,{'integrate_center' :  indx[k], 'integrate_width' : 10})
        if len(iterativeopt_workspace['proc'].values) > 1:
            intgrl_array.append(sum(abs(iterativeopt_workspace['proc'].values)))
        else:
            intgrl_array.append(abs(iterativeopt_workspace['proc'].values[0]))
    cent = np.argmax(intgrl_array)
    
    return indx[cent]

total_folders = []
total_folders.append(folder_p0)
for e in folders_Enhancements:
    total_folders.append(e)
for t in folders_T1s:
    total_folders.append(t)
total_folders.append(folder_T10)

T1 = []
T1_stdd = []
E = []
for f in range(0, len(total_folders)):

    data = dnp.dnpImport.topspin.import_topspin(directory  + os.sep,total_folders[f])
    workspace = dnp.create_workspace('raw',data)
    workspace.copy('raw','proc')

    dnp.dnpNMR.remove_offset(workspace,{})
    dnp.dnpNMR.window(workspace,{'linewidth' : proc_params['window_linewidth']})
    dnp.dnpNMR.fourier_transform(workspace,{'zero_fill_factor' : proc_params['zero_fill_factor']})
    
    if workspace['proc'].ndim == 2:
        workspace = dnp.dnpNMR.autophase(workspace,{})
        workspace = dnp.dnpNMR.align(workspace, {})
    else:
        phase = workupPhaseOpt(workspace)
        workspace['proc'] *= np.exp(-1j * phase)
     
    int_params['integrate_center'] = optCenter(workspace)
    workspace = dnp.dnpNMR.integrate(workspace,{'integrate_center' :  int_params['integrate_center'], 'integrate_width' : int_params['integrate_width']})

    if  total_folders[f] == folder_p0:
        p0 = np.real(workspace['proc'].values[0])
        print('Done with p0 folder...')
    elif  total_folders[f] == folder_T10:
        workspace = dnp.dnpFit.t1Fit(workspace)
        T10 = workspace['fit'].attrs['t1']
        T10_stdd = workspace['fit'].attrs['t1_stdd']
        print('Done with T1(0) folder...')
    elif total_folders[f] in folders_T1s:
        workspace = dnp.dnpFit.t1Fit(workspace)
        T1.append(workspace['fit'].attrs['t1'])
        T1_stdd.append(workspace['fit'].attrs['t1_stdd'])
        print('Done with T1(p) folder ' + str(total_folders[f]) + '...')
    elif total_folders[f] in folders_Enhancements:
        E.append(np.real(workspace['proc'].values[0]) / p0)
        print('Done with Enhancement(p) folder ' + str(total_folders[f]) + '...')

E = np.array([Enhancement_powers, E])
E = np.transpose(E)
E = E[E[:,0].argsort()]

T1 = np.array([T1_powers, T1])
T1 = np.transpose(T1)
T1 = T1[T1[:,0].argsort()]

Enhancement_powers = E[:,0]
Enhancements = E[:,1]

T1_powers = T1[:,0]
T1s = T1[:,1]

print('Passing to dnpHydration...')

hydration = {
'E' : np.array(Enhancements),
'E_power' : np.array(Enhancement_powers),
'T1' : np.array(T1s),
'T1_power' : np.array(T1_powers)
}

hydration.update({
'T10': T10,
'T100': T100,
'spin_C': spin_C,
'field': field,
'smax_model': smax_model,
't1_interp_method': t1_interp_method
})

hydration_workspace = dnp.create_workspace()
hydration_workspace.add('hydration_inputs', hydration)
                     
hydration_results = dnp.dnpHydration.hydration(hydration_workspace)

hydration_workspace.add('hydration_results', hydration_results)
hydration_workspace['hydration_results'].update({'T1p_stdd': T1_stdd, 'T10_stdd': T10_stdd})

svpthnm = directory + save_name
if os.path.isfile(directory + save_name + '.h5'):
    svpthnm = svpthnm + '_COPY'
    
dnp.dnpImport.h5.saveh5(hydration_workspace, svpthnm + '.h5')

print('Done with dnpHydration.')
print('T1(0) = ' + str(T10) + ', stdd = ' + str(T10_stdd))
print('krho = ' + str(hydration_results['krho']))
print('ksigma = ' + str(hydration_results['ksigma']) + ', stdd = ' + str(hydration_results['ksigma_stdd']))
print('coupling factor = ' + str(hydration_results['coupling_factor']))
print('tcorr = ' + str(hydration_results['tcorr']))
