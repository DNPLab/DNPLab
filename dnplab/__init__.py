# sub-modules should be imported first
from . import (
    core,
    dnpFit,
    dnpHydration,
    dnpIO,
    dnpNMR,
    dnpResults,
    dnpTools,
    hydrationGUI,
)

# before raising the functions below in the sub-modules to the module
from .dnpData import create_workspace, dnpdata, dnpdata_collection
from .dnpImport import load
from .version import __version__
