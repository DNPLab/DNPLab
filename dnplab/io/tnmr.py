from .. import DNPData
import numpy as _np
import struct

TMAG_ATTRS ={
    # format:
    # attrs: ('dtype', start_index, end_index) 
    # Number of points and scans in all dimensions:
    'npts': ('<4i', 0, 16),
    'acutual_npts': ('<4i', 16, 32),
    'acq_points': ('<i', 32, 36),
    'npts_start': ('<4i', 36, 52),
    'scans': ('<i', 52, 56),
    'actual_scans': ('<i', 56, 60),
    'dummy_scans': ('<i', 60, 64),
    'repeat_times': ('<i', 64, 68),
    'sadimension': ('<i', 68, 72),
    'samode': ('<i', 72, 76),

    # Field and frequencies:
    'magnet_field': ('<d', 76, 84),
    'ob_freq': ('<4d', 84, 116),
    'base_freq': ('<4d', 116, 148),
    'offset_freq': ('<4d', 148, 180),
    'ref_freq': ('<d', 180, 188),
    # 'nmr_freq': ('<d', 188, 196),
    'obs_channel': ('<h', 196, 198),
    # space from 198 to 240

    # Spectral width, dwell and filter:
    'sw': ('<4d', 240, 272),
    'dwell': ('<4d', 272, 304),
    'filter': ('<d', 304, 312),
    'experiment_time': ('<d', 312, 320),
    'acq_time': ('<d', 320, 328),
    'last_delay': ('<d', 328, 336),
    'spectrum_direction': ('<h', 336, 338),
    'hardware_sideband': ('<h', 338, 340),
    'taps': ('<h', 340, 342),
    'type': ('<h', 342, 344),
    'bDigRec': ('<i', 344, 348),
    'nDigitalCenter': ('<i', 348, 352),
    # space from 352 to 368

    # Hardware settings: 
    'transmitter_gain': ('<h', 368, 370),
    'receiver_gain': ('<h', 370, 372),
    'NumberOfReceivers': ('<h', 372, 374),
    'RG2': ('<h', 374, 376),
    'receiver_phase': ('<d', 376, 384),
    # space from 384 to 388

    # Spinning speed information:
    'set_spin_rate': ('<H', 388, 390),
    'actual_spin_rate': ('<H', 390, 392),

    # Lock information:
    'lock_field': ('<h', 392, 394),
    'lock_power': ('<h', 394, 396),
    'lock_gain': ('<h', 396, 398),
    'lock_phase': ('<h', 398, 400),
    'lock_freq_mhz': ('<d', 400, 408),
    'lock_ppm': ('<d', 408, 416),
    'H2O_freq_ref': ('<d', 416, 424),
    # space from 424 to 440

    # VT information:
    'set_temperature': ('<d', 440, 448),
    'actual_temperature': ('<d', 448, 456),

    # Shim information:
    'shim_units': ('<d', 456, 464),
    'shims': ('<36h', 464, 536),
    'shim_FWHM': ('<d', 536, 544),

    # Bruker specific information:
    'HH_dcpl_attn': ('<h', 544, 546),
    'DF_DN': ('<h', 546, 548),
    'F1_tran_mode': ('<7h', 548, 562),
    'dec_BW': ('<h', 562, 564),
    'grd_orientation': ('<4c', 564, 568),
    'LatchLP': ('<i', 568, 572),
    'grd_Theta': ('<d', 572, 580),
    'grd_Phi': ('<d', 580, 588),
    # space from 588 to 852

    # Time variables:
    'start_time': ('<i', 852, 856),
    'finish_time': ('<i', 856, 860),
    'elapsed_time': ('<i', 860, 864),

    # Text variables:
    'date': ('<32c', 864, 896),
    'nucleus': ('<16c', 896, 912),
    'nucleus_2D': ('<16c', 912, 928),
    'nucleus_3D': ('<16c', 928, 944),
    'nucleus_4D': ('<16c', 944, 960),
    'sequence': ('<32c', 960, 992),
    'lock_solvent': ('<16c', 992, 1008),
    'lock_nucleus': ('<16c', 1008, 1024),
}

def import_tnmr(path, squeeze=True):
    """Import tnmr data and return DNPData object

    Args:
        path (str):         Path to .jdf file
        squeeze (bool):     Automatically remove length 1 dimensions

    Returns:
        dnpdata (object):   DNPData object containing tnmr data
    """

    values, dims, coords, attrs = import_tnmr_data(path)

    tnmr_data = DNPData(values, dims, coords, attrs)

    if squeeze:
        tnmr_data.squeeze()

    return tnmr_data


def import_tnmr_data(path):
    """Import spectrum or spectra of tnmr data

    Args:
        path (str):         Path to .tnt file

    Returns:
        data (ndarray):     Spectrum or spectra if >1D
        abscissa (list):    Coordinates of axes
        dims (list):        Axes names
    """

    with open(path, "rb") as f:
        raw = f.read()

    tmag_attrs = {}
    tmag_attrs["version"] = str(raw[0:8].decode("utf-8")).replace('\x00', '')
    attrs["experiment_type"] = "nmr_spectrum"

    len_tecmag_struct = int.from_bytes(raw[16:20], "little")
    tecmag_struct = raw[20 : 20 + len_tecmag_struct]

    offset = 28 + len_tecmag_struct
    len_data = raw[offset : offset + 4]
    offset += 4

    len_data = int.from_bytes(len_data, "little")

    data = raw[offset : offset + (len_data)]

    data = _np.frombuffer(data, dtype="<f")

    data = data[0::2] + 1j * data[1::2]

    # Parse tecmag struct
    for key, val in TMAG_ATTRS.items():
        if 'c' not in val[0]:
            tmag_attrs[key] = struct.unpack(val[0], tecmag_struct[val[1]:val[2]])
        else:
            tmag_attrs[key] = str(tecmag_struct[val[1]:val[2]].decode("utf-8")).replace('\x00', '')

    # points in x0, x1, x2, x3
    npts = tmag_attrs['npts']

    # Array of points in each dimension (x0, x1, x2, x3)
    tmag_attrs['nmr_frequency'] = tmag_attrs['ob_freq'][0] * 1e6 # nmr_frequency is equal to ob_freq in MHz

    dwell = tmag_attrs['dwell']

    data = data.reshape(npts, order="F")

    coords = []

    for ix, pts in enumerate(npts):
        coord = _np.r_[0:pts] * dwell[ix]
        coords.append(coord)

    dims = ["t2", "t1", "t3", "t4"]  # t2 dim is first

    return data, dims, coords, attrs
