import sys
sys.path.append('..')
import dnplab

exp_num = 40
path = '../data/prospa/toluene_10mM_Tempone/%i/'

data = dnplab.dnpImport.prospa.import_prospa(path%exp_num)

ws = dnplab.create_workspace('raw', data)
ws.copy('raw')

ws = dnplab.dnpNMR.window(ws, {})
ws = dnplab.dnpNMR.fourier_transform(ws, {})
ws = dnplab.dnpNMR.autophase(ws, {})
ws.copy('proc', 'ft')
ws = dnplab.dnpNMR.integrate(ws, {})

dnplab.dnpResults.plot(ws['ft'].real, label = 'toluene')
dnplab.dnpResults.legend()
#dnplab.plot(ws['proc'].imag)
dnplab.dnpResults.xlim(20,-20)
dnplab.dnpResults.show()
