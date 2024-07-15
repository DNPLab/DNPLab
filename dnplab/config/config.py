"""
global config
"""

import configparser
from pathlib import Path
import warnings

import logging

logger = logging.getLogger(__name__)


def _escape_split(s, delim=",", escape="\\"):
    tokens = []
    previous_escape = False
    subtoken = ""
    for k in range(len(s)):
        if s[k] == delim and (not previous_escape):
            if len(subtoken) > 0:
                tokens.append(subtoken)
            subtoken = ""  # reset subtoken
        else:
            # ESCAPE DELIM  -> DELIM
            if previous_escape and s[k] != escape and s[k] == delim:
                subtoken = subtoken[:-1] + s[k]
            else:
                # add to subtoken
                subtoken += s[k]
            # set previous_escape flag:
            # If for current char is escape (True and s[k]=='\\') and previous_escape is False -> set it to True
            # If for current char is escape (True and s[k]=='\\') and previous_escape is True -> case of '\\\\' -> escaping an escape character -> set it back to False
            # if current char is no escape character -> set it to false
            previous_escape = (not previous_escape) and (s[k] == escape)
    if len(subtoken) > 0:
        tokens.append(subtoken)
    return tokens


def _kwarg_converter(s: str):
    tokens = _escape_split(s, ",", escape="\\")
    args = []
    kwargs = {}
    for k in tokens:
        subtokens = _escape_split(k, "=", escape="\\")
        if len(subtokens) == 1:
            args.append(subtokens[0])
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
