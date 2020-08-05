import sys
sys.path.append('..')
import numpy as np

import dnpLab as dnp

path = '../data/prospa/toluene_10mM_Tempone/'
exp_list = [str(x) for x in range(1,43)]
#exp_list = ['1','2','3','4']
power_list = 10.**(np.r_[0:40:len(exp_list)*1j] / 10) / 1000.
ws = dnp.dnpImport.prospa.import_prospa_dir(path, exp_list)

for ix, exp_num in enumerate(exp_list):
    tmp = ws[exp_num]
    tmp.new_dim('power', power_list[ix])
    if ix == 0:
        data = tmp
    else:
        data.concatenate(tmp, 'power')

ws.add('raw', data)
ws.copy('raw', 'proc')

ws = dnp.dnpNMR.window(ws, {})
ws = dnp.dnpNMR.fourier_transform(ws, {})
ws['proc'] = ws['proc']['t2',(-100,100)]
ws = dnp.dnpNMR.align(ws, {})
ws.copy('proc', 'ft')
ws = dnp.dnpNMR.integrate(ws, {'integrate_width': 40})

dnp.plot(ws['ft'])
dnp.xlim(20,-20)
dnp.figure()
dnp.plot(ws['proc'], 'o')
dnp.show()

