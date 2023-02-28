import pathlib
import unittest
import dnplab as dnp
import os

import numpy as np

class helper_functions_processing(unittest.TestCase):
    def setUp(self):
        self.testdata = os.path.join(".", "data", "csv")
        p = pathlib.Path(self.testdata)
        self.data = dnp.io.load_csv.load_csv(
            p.joinpath("csv_example.csv"),
            skiprows=1,
            maxrows=1000,
            tcol=0,
            real=1,
            imag=3,
        )

    def test_000_funcionality_signal_to_noise(self):
        """
        check only whether the function raises no error with DNPData input, not whether rsults are useful
        alot of simple tests lumped together

        Missing: check whether signal and noise are scalar values
        """
        f= dnp.processing.signal_to_noise

        self.assertRaises(ValueError,f,self.data,(300,400),(500,600) )
        data=dnp.fourier_transform(self.data)

        try:
            snr=f(data,(300,400),(500,600) )
        except ValueError as e:
            self.fail('signal_to_noise reported ValueError {0}'.format(e))
        self.assertTrue(not np.isnan(snr))

        snr,signal,noise=f(data,(300,400),(500,600),fullOutput=True,detrend=False)

        self.assertTrue( not np.isnan(snr))
        self.assertTrue( not np.isnan(signal))
        self.assertTrue( not np.isnan(noise))



if __name__ == "__main__":
    unittest.main()
    pass
