import unittest
from numpy.testing import assert_array_equal
import dnplab as dnp
import numpy as np

class expoinential_window_tester(unittest.TestCase):
    def setUp(self):
        lw = 1
        x = np.r_[0:10:1024j]
        y = dnp.window.exponential(x, lw)
        self.data = dnp.DNPData(y, ['t2'], [x])
        
        dt = x[1] - x[0]
        self.f = np.r_[-1/(2*dt):1/(2*dt):1024j]
        self.lorizen_distribution = 1 / dnp.pi * (1 / ((self.f)**2 + 0.5*lw**2))
        self.lorizen_distribution /= np.max(self.lorizen_distribution)
        print(self.lorizen_distribution)
        
    def test_exponential_window(self):
        self.data = dnp.fourier_transform(self.data, convert_to_ppm=False)
        self.data /= np.max(self.data)
        print(self.data.values)
        dnp.plt.figure()
        dnp.plot(self.data)
        dnp.plt.plot(self.f, self.lorizen_distribution)
        dnp.plt.show()

if __name__ == "__main__":
    unittest.main()
