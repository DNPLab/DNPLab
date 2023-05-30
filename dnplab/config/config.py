"""
global config
"""
import configparser
from pathlib import Path
import warnings


def _kwarg_converter(s: str):
    tokens = s.strip("[").strip("]").split(",")
    args = []
    kwargs = {}
    for k in tokens:
        subtokens = k.split("=")
        if len(subtokens) == 1:
            args.append(k[0])
        else:
            kwargs[subtokens[0].strip()] = subtokens[1].strip()
    return args, kwargs


def _get_dnp_config(configname="dnplab.cfg"):
    config = configparser.ConfigParser(
        converters={
            "list": lambda x: list(x.strip("[").strip("]").split(",")),
            "args_kwargs": _kwarg_converter,
        }
    )

    # define three possible locations:
    dnplab_current_config = Path.cwd() / configname
    dnplab_home_config = Path.home() / configname

    dnplab_cfg_folder = str(
        Path(__file__).parent
    )  # / configname #.with_name("config"))
    dnplab_global_config = Path(dnplab_cfg_folder) / configname

    config_read_list = [dnplab_global_config, dnplab_home_config, dnplab_current_config]

    # user defined takes precedence
    config.read(config_read_list)
    return config


DNPLAB_CONFIG = _get_dnp_config()
