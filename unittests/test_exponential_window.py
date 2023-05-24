import unittest
from numpy.testing import assert_array_equal
import dnplab as dnp
import numpy as np

class exponential_window_tester(unittest.TestCase):
    def setUp(self):
        lw = 1
        x = np.r_[0:10:1024j]
        y = dnp.window.exponential(x, lw)
        y = y + 0j
        
        self.data = dnp.DNPData(y, ['t2'], [x])
        self.data_ft = dnp.fourier_transform(self.data.copy(), convert_to_ppm=False, zero_fill_factor = 1)
        self.data_ft /= np.max(self.data_ft)
        self.f = self.data_ft.coords['f2']
        self.lorentzian_distribution = 1 / dnp.pi * (1 / ((self.f)**2 + 0.5*lw**2))
        self.lorentzian_distribution /= np.max(self.lorentzian_distribution)
        
    def test_exponential_window(self):
        dnp.plt.figure()
        dnp.plot(self.data, label = 'window')

        dnp.plt.figure()
        dnp.plot(self.data_ft, label = 'FFT window')
        dnp.plt.plot(self.f, self.lorentzian_distribution, label = 'Lorentzian')
        dnp.plt.legend()
        dnp.plt.show()

if __name__ == "__main__":
    unittest.main()