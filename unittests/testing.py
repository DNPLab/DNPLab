import numpy as np
import dnplab as dnp


def get_gauss_3d(std_noise=0.0):
    x = np.r_[0:100]
    y = np.r_[0:100]
    z = np.r_[0:100]

    noise = std_noise * np.random.randn(len(x), len(y), len(z))
    gauss = np.exp(-1.0 * (x - 50) ** 2.0 / (10.0 ** 2))
    gauss_3d = (
        gauss.reshape(-1, 1, 1) * gauss.reshape(1, -1, 1) * gauss.reshape(1, 1, -1)
    )
    gauss_3d += noise
    return x, y, z, gauss_3d


def test3d(std_noise=0.0):
    x, y, z, gauss_3d = get_gauss_3d(std_noise)
    test_data = dnp.DNPData(gauss_3d,  ["x", "y", "z"], [x, y, z])
    return test_data
