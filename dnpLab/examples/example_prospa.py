import sys
sys.path.append('../..')
import dnpLab as dnp


path = '../../data/kea/toluene_10mM_Tempone/30/'

data = dnp.dnpImport.kea.importKea(path)

ws = dnp.create_workspace('raw', data)
ws.copy('raw')

ws = dnp.dnpNMR.window(ws, {})
ws = dnp.dnpNMR.fourier_transform(ws, {})
ws = dnp.dnpNMR.autophase(ws, {})

dnp.plot(ws['proc'].real)
dnp.plot(ws['proc'].imag)
dnp.show()
