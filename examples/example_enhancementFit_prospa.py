import sys
sys.path.append('..')

import numpy as np
import dnpLab as dnp

path = '../data/prospa/toluene_10mM_Tempone/'
exp_list = [str(x) for x in range(1,43)]
#exp_list = ['1','2','3','4']
power_list = 10.**(np.r_[0:40:len(exp_list)*1j] / 10) / 1000.

for ix, exp_num in enumerate(exp_list):
    tmp = dnp.dnpImport.prospa.import_prospa(path + '%s'%exp_num)
    tmp.new_dim('power', power_list[ix])
    if ix == 0:
        data = tmp
    else:
        data.concatenate(tmp, 'power')

ws = dnp.create_workspace()
ws.add('raw', data)
ws.copy('raw', 'proc')

ws = dnp.dnpNMR.window(ws, {})
ws = dnp.dnpNMR.fourier_transform(ws, {})
ws['proc'] = ws['proc']['t2',(-100,100)]
ws = dnp.dnpNMR.align(ws, {})
ws.copy('proc', 'ft')
ws = dnp.dnpNMR.integrate(ws, {'integrate_width': 40})

ws = dnp.dnpFit.enhancementFit(ws)

dnp.dnpResults.plot(ws['ft'])
dnp.dnpResults.xlim(20,-20)
dnp.dnpResults.figure()
dnp.dnpResults.plot(ws['proc'], 'o')
dnp.dnpResults.plot(ws['fit'])
dnp.dnpResults.show()

