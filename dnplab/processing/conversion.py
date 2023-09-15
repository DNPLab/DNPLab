import numpy as _np
from ..core.data import DNPData


def convert_power(data, mode="dBm2W"):
    """Convert power in dBm to power in W and vice versa

    Convert power in dBm to power in W and vice versa

    If power_unit exists in the DNPData object, the function will automatically convert the power, regardless the mode.

    If power_unit doesn't exist, the function will be performed in the default mode.

    If data is not DNPLab object, the function will be performed based on the mode.

    Args:
        data (DNPData, array, list, float or int)
        mode (str): By default, the function convert power in dBm to power in W

    Returns:
        out (DNPData, array, list, float or int)

    """
    if isinstance(data, DNPData):
        out = data.copy()
        dims = out.dims
        for dim in dims:
            if dim.lower() == "power" or dim.lower() == "powers":
                break
            elif dim == dims[-1]:
                raise Warning("Power is not in dim")
            else:
                continue

        if "power_unit" in out.dnplab_attrs.keys():
            if out.dnplab_attrs["power_unit"] == "dBm":
                mode = "dBm2W"
            elif out.dnplab_attrs["power_unit"] == "W":
                mode = "W2dBm"
            else:
                raise Warning("Power unit in dnplab_attrs is invalid")

        if mode.lower() == "dbm2w":
            power_unit = "W"
            f = dBm2w
        elif mode.lower() == "w2dbm":
            power_unit = "dBm"
            f = w2dBm
        else:
            raise Warning("Mode is not acceptable")

        out.coords[dim] = f(out.coords[dim])
        proc_parameters = {
            "mode": mode,
        }
        out.dnplab_attrs["power_unit"] = power_unit
        proc_attr_name = "convert_power"
        out.add_proc_attrs(proc_attr_name, proc_parameters)

        return out

    else:
        if mode.lower() == "dbm2w":
            f = dBm2w
        elif mode.lower() == "w2dbm":
            f = w2dBm
        else:
            raise Warning("Mode is not acceptable")
        return f(data)


def dBm2w(power_in_dBm):
    """Convert power in dBm to power in W

    Convert a microwave power given in dBm to W. Note that for values <= -190 dBm the output power in W is set to 0.

    Args:
        power_in_dBm (array, list, float or int): Power in (dBm)

    Returns:
        float (array, list, float or int): Power in (W)

    """

    if isinstance(power_in_dBm, (_np.ndarray, list)):
        power_in_W = power_in_dBm.copy()
        for index in range(len(power_in_dBm)):
            power_in_W[index] = dBm2w(power_in_dBm[index])
        return power_in_W

    else:
        power_in_W = 10.0 ** (power_in_dBm / 10.0) / 1000.0
        # Set values below 1 pW to 0 W
        power_in_W = 0 if power_in_W <= 1e-22 else power_in_W

    return power_in_W


def w2dBm(power_in_W):
    """Convert power in W to power in dbM

    Convert a microwave power given in W to dBm

    Args:
        power_in_W (array, list, float or int):   Power in (W)

    Returns:
        float (array, list, float or int): Power in (W)

    """
    if isinstance(power_in_W, (_np.ndarray, list)):
        power_in_dBm = power_in_W.copy()
        for index in range(len(power_in_W)):
            power_in_dBm[index] = w2dBm(power_in_W[index])
        return power_in_dBm

    else:
        power_in_dBm = 10.0 * _np.log10(1000 * power_in_W)

    return power_in_dBm


def calc_tp90(c, P, Q=1, alpha=0, verbose=False):
    """Calculate 90 degree pulse length

    Calculate 90 degree pulse length (tp90) from probe conversion factor, and applied RF power. Optionally, the quality factor and attenuation can be given as input arguments. The function returns the pulse length in (ns). A formatted output can be generated when setting the verbose flag to True.

    Args:
        c (float):          Probe conversion factor (G/sqrt(W))
        P (float):          Input RF power
        Q (float):          Optionally, probe quality factor. Default value is 1
        alpha (float):      Optionally, attenuation (dB)
        verbose (boolean):  Optionally, return results in formatted output

    Returns:
        tp90 (float):       90 degree pulse length (ns)

    .. math::

    """

    power_at_probe = P / (10 ** (alpha / 10))
    b1_g = c * _np.sqrt(power_at_probe * Q)
    b1_mhz = b1_g * 2.804 * 1e-6

    tp90 = 1 / b1_mhz / 4 * 1e-12

    if verbose == True:
        print(" ")
        print("*** Input Parameters ***")
        print("Conversion Factor c (G/sqrt(W)): ", c)
        print("Input RF Power P (W):            ", P)
        print("Quality Factor Q:                ", Q)
        print("Attenuation alpha (dB):          ", alpha)
        print(" ")
        print("*** Results ***")
        print("RF power at probe (W):           ", power_at_probe)
        print("B1 Field Strength (G):           ", b1_g)
        print("B1 Field Strength (MHz):         ", b1_mhz * 1e6)
        print("tp90 (ns):                       ", tp90 * 1e9)

    return tp90


def calc_conversion_factor(tp90, P, Q=1, alpha=0, verbose=False):
    """Calculate probe conversion factor

    Calculate the probe microwave conversion factor from the 90 degree pulse length (tp90), and applied RF power. Optionally, the quality factor and attenuation can be given as input arguments. The function returns the conversion factor in (G/sqrt(W)). A formatted output can be generated when setting the verbose flag to True.

    Args:
        tp90 (float):       90 degree pulse length (ns)
        P (float):          Input RF power
        Q (float):          Optionally, probe quality factor. Default value is 1
        alpha (float):      Optionally, attenuation (dB)
        verbose (boolean):  Optionally, return results in formatted output

    Returns:
        c (float):          Probe conversion factor (G/sqrt(W))

    .. math::

    """

    power_at_probe = P / (10 ** (alpha / 10))

    b1_mhz = 1 / tp90 / 4 * 1e-12
    b1_g = b1_mhz / 2.804e-6

    c = b1_g / _np.sqrt(power_at_probe * Q)

    if verbose == True:
        print(" ")
        print("*** Input Parameters ***")
        print("tp90 (ns):                       ", tp90 * 1e9)
        print("Input Power P (W):               ", P)
        print("Quality Factor Q:                ", Q)
        print("Attenuation alpha (dB):          ", alpha)
        print(" ")
        print("*** Results ***")
        print("Power at probe (W):              ", power_at_probe)
        print("B1 Field Strength (G):           ", b1_g)
        print("B1 Field Strength (MHz):         ", b1_mhz * 1e6)
        print("Conversion Factor c (G/sqrt(W)): ", c)
