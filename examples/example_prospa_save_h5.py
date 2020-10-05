import sys

sys.path.append("..")
import numpy as np

import dnplab

path = "../data/prospa/toluene_10mM_Tempone/"
exp_list = [str(x) for x in range(1, 43)]
# exp_list = ['1','2','3','4']
power_list = 10.0 ** (np.r_[0 : 40 : len(exp_list) * 1j] / 10) / 1000.0
tmp_ws = dnplab.dnpImport.prospa.import_prospa_dir(path, exp_list)
ws = dnplab.create_workspace()

for ix, exp_num in enumerate(exp_list):
    tmp = tmp_ws[exp_num]
    tmp.new_dim("power", power_list[ix])
    if ix == 0:
        data = tmp
    else:
        data.concatenate(tmp, "power")

ws.add("raw", data)
ws.copy("raw", "proc")

ws = dnplab.dnpNMR.window(ws, {})
ws = dnplab.dnpNMR.fourier_transform(ws, {})
ws["proc"] = ws["proc"]["t2", (-100, 100)]
ws = dnplab.dnpNMR.align(ws, {})
ws.copy("proc", "ft")
ws = dnplab.dnpNMR.integrate(ws, {"integrate_width": 40})

dnplab.dnpResults.plot(ws["ft"])
dnplab.dnpResults.xlim(20, -20)
dnplab.dnpResults.figure()
dnplab.dnpResults.plot(ws["proc"], "o")

dnplab.dnpImport.h5.saveh5(ws, "test.h5")

dnplab.dnpResults.show()
