def attrs4dnplab(exp_attrs):
    """Convert experiment attributes to dnplab attributes

    Args:
        exp_attrs (dict): Dictionary of prospa experiment acqusition parameters

    Returns:
        dnplab_attrs (dict): Dictionary of parameters used in dnplab
    """
    dnplab_attrs = {}
    dnplab_attrs["experiment_type"] = "nmr_spectrum"
    dnplab_attrs["spectrometer_format"] = exp_attrs["spectrometer_format"]
    dnplab_attrs["spectrometer_frequency"] = exp_attrs["nmr_frequency"] # Hz

    if exp_attrs["spectrometer_format"] == "prospa":

        if '1Pulse' in exp_attrs["experiment"]:
            dnplab_attrs["experiment"] = "1D"
    
        elif "T1-IR-FID" in exp_attrs["experiment"]:
            dnplab_attrs["experiment"] = "2D IR"
            dnplab_attrs["minimum_delay"] = exp_attrs["minDelay"] * 1e-3 # s
            dnplab_attrs["maximum_delay"] = exp_attrs["maxDelay"] * 1e-3 # s
        
        elif "jres2D" in exp_attrs["experiment"]:
            dnplab_attrs["experiment"] = "2D JRES"
            dnplab_attrs["number_of_steps"] = exp_attrs["nrSteps"]
            dnplab_attrs["inter_pulse_delay"] = exp_attrs["interPulseDelay"]
            dnplab_attrs["increment"] = exp_attrs["increment"]

        else:
            dnplab_attrs["experiment"] = "Not Defined"
    
        dnplab_attrs["number_of_scans"] = exp_attrs["nrScans"]
        dnplab_attrs["90_pulse_length"] = exp_attrs["90Amplitude"]
        dnplab_attrs["receiver_gain"] = exp_attrs["rxGain"]


    return dnplab_attrs