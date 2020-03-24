from biolink.dblink import KEGGLinker, GeneNameLinker, DrugBankLinker, SiderLinker, UniprotLinker
from os.path import realpath, dirname

file_dp = dirname(realpath(__file__))
__LIB_ROOT__DIR__ = dirname(file_dp)

__name__ = "biolink"
__description__ = "A library for linking entities of biological knowledge bases."
__version__ = "0.0.1"
__all__ = ["GeneNameLinker", "KEGGLinker", 'DrugBankLinker', 'SiderLinker', 'UniprotLinker', "__version__", "__name__", "__description__"]
