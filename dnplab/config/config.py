"""
global config
"""
import configparser
from pathlib import Path
import warnings


def _get_dnp_config(configname="dnplab_cfg.cfg"):
    config = configparser.ConfigParser(
        converters={
            "list": lambda x: list(
                x.strip("(").strip("[").strip("]").strip(")").split(",")
            )
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


if not "DNPLAB_CONFIG" in locals():
    DNPLAB_CONFIG = _get_dnp_config()
