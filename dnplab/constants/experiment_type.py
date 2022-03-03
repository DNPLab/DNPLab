"""Dictionary to store experiment types"""

experiment_type_list =[
    "nmr_pectrum",
    "epr_pectrum",
    "integrals",
    "enhancements",
    "inversion_recovery"
]





experiment_type_dict = {}

experiment_type_dict["nmr"] = ["spectrum", "fid"]
experiment_type_dict["epr"] = ["cw_spectrum", "echo_detected_spectrum", "deer_time_trace", "eseem_time_trace"]
experiment_type_dict["processing"] = ["integrals", "enhancements"]


