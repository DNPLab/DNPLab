def attrs4dnplab(exp_attrs):
    """Convert experiment attributes to dnplab attributes

    Args:
        exp_attrs (dict): Dictionary of prospa experiment acqusition parameters

    Returns:
        dnplab_attrs (dict): Dictionary of parameters used in dnplab
    """
    dnplab_attrs = {}
    dnplab_attrs["experiment_type"] = exp_attrs["experiment_type"]
    dnplab_attrs["spectrometer_format"] = exp_attrs["spectrometer_format"]

    if exp_attrs["spectrometer_format"] == "prospa":
        dnplab_attrs["spectrometer_frequency"] = exp_attrs["nmr_frequency"]  # Hz

        if "1Pulse" in exp_attrs["experiment"]:
            dnplab_attrs["experiment"] = "1D"

        elif "T1-IR-FID" in exp_attrs["experiment"]:
            dnplab_attrs["experiment"] = "2D IR"
            dnplab_attrs["minimum_delay"] = exp_attrs["minDelay"] * 1e-3  # s
            dnplab_attrs["maximum_delay"] = exp_attrs["maxDelay"] * 1e-3  # s

        elif "jres2D" in exp_attrs["experiment"]:
            dnplab_attrs["experiment"] = "2D JRES"
            dnplab_attrs["number_of_steps"] = exp_attrs["nrSteps"]
            dnplab_attrs["inter_pulse_delay"] = exp_attrs["interPulseDelay"]
            dnplab_attrs["increment"] = exp_attrs["increment"]

        else:
            dnplab_attrs["experiment"] = "Not Defined"
        dnplab_attrs["spectrometer_frequency"] = exp_attrs["nmr_frequency"]  # Hz
        dnplab_attrs["90_pulse_length"] = exp_attrs["90Amplitude"]
        dnplab_attrs["receiver_gain"] = exp_attrs["rxGain"]

    elif exp_attrs["spectrometer_format"] == "topspin":
        dnplab_attrs["experiment"] = exp_attrs["experiment"]

    elif exp_attrs["spectrometer_format"] == "xepr":
        dnplab_attrs["microwave_frequency"] = exp_attrs["frequency"] * 1e9
        dnplab_attrs["experiment"] = exp_attrs["experiment"]

        if exp_attrs["x_unit"] == "G":
            dnplab_attrs["center_field"] = exp_attrs["center_field"] * 1e5  # T
        elif exp_attrs["x_unit"] == "T":
            dnplab_attrs["center_field"] = exp_attrs["center_field"]  # T
        if dnplab_attrs["experiment"] != "2D ESE":
            dnplab_attrs["power"] = exp_attrs["power"]  # W
            if dnplab_attrs["experiment"] in ["1D", "2D"]:
                dnplab_attrs["attenuation"] = exp_attrs["attenuation"]
                dnplab_attrs["conversion_time"] = (
                    exp_attrs["conversion_time"] * 1e-3
                )  # s
                dnplab_attrs["modulation_amplitude"] = exp_attrs["modulation_amplitude"]
                dnplab_attrs["modulation_frequency"] = (
                    exp_attrs["modulation_frequency"] * 1e3
                )  # Hz
                dnplab_attrs["time_constant"] = exp_attrs["time_constant"] * 1e-3  # s
        else:
            dnplab_attrs["pulse_attenuation"] = exp_attrs["pulse_attenuation"]
    dnplab_attrs["number_of_scans"] = exp_attrs["nrScans"]

    return dnplab_attrs
