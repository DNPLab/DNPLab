"""DNPLab - Bringing the Power of Python to DNP-NMR Spectroscopy"""

from .core.data import DNPData
from .core.ufunc import *
from .core.util import *

from .constants import *
from .fitting import *
from .math import *

from .io import *
from .io.save import save
from .io.load import load

from .analysis.relaxation_fit import *
from .analysis.hydration import hydration
from .analysis.simulate_enhancement_profiles import *

from .processing import *
from .widgets import *
from .plotting import *
from .reporting import *
from .version import __version__


"""
v2 test, atm keine kommentare
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

    dnplab_cfg_folder = str(Path(__file__).with_name("config"))
    dnplab_global_config = Path(dnplab_cfg_folder) / configname

    config_read_list = [dnplab_global_config, dnplab_home_config, dnplab_current_config]

    # user defined takes precedence
    config.read(config_read_list)
    return config


DNPLAB_CONFIG = _get_dnp_config()
