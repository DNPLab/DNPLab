import sys
sys.path.append('../..')
import dnpLab as dnp

exp_num = 40
path = '../../data/prospa/toluene_10mM_Tempone/%i/'

data = dnp.dnpImport.prospa.import_prospa(path%exp_num)

ws = dnp.create_workspace('raw', data)
ws.copy('raw')

ws = dnp.dnpNMR.window(ws, {})
ws = dnp.dnpNMR.fourier_transform(ws, {})
ws = dnp.dnpNMR.autophase(ws, {})
ws.copy('proc', 'ft')
ws = dnp.dnpNMR.integrate(ws, {})

dnp.plot(ws['ft'].real, label = 'toluene')
dnp.legend()
#dnp.plot(ws['proc'].imag)
dnp.xlim(20,-20)
dnp.show()
