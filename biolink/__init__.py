from biolink.dblink import *
from os.path import realpath, dirname

file_dp = dirname(realpath(__file__))
__LIB_ROOT__DIR__ = dirname(file_dp)

__name__ = "biolink"
__description__ = "A library for linking entities of biological knowledge bases."
__version__ = "0.0.2"
__all__ = [
    "GeneNameLinker", "KEGGLinker", 'DrugBankLinker', 'SiderLinker',
    'UniprotLinker', "BiogridLinker", "ChemblLinker", "DipLinker",
    "EmblLinker", "EmblCDSLinker", "EnsemblLinker", "GeneDBLinker",
    "GILinker", "HGNCLinker", "HPALinker", "MimLinker",
    "MintLinker", "PharmgkbLinker", "ProteomicsDBLinker", "RefSeqLinker",
    "StringLinker", "UniparcLinker", "UniprotKBLinker", "Uniref100Linker",
    "ThreeDMetLinker", "ChemspiderLinker", "HMDBLinker", "GuidePharmaLinker",
    "PubchemLinker", "ChebiLinker", "DPDLinker", "BindingDBLinker",
    "HSDBLinker", "IupharLinker", "KnapsackLinker", "LingandBoxLinker",
    "MassBankLinker", "NikkajiLinker", "PDBLinker", "TTDLinker",
    "WikipediaLinker", "CellosaurusLinker",
    "__version__", "__name__", "__description__"
]
