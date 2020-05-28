from biodblinker.dblink import *
from biodblinker.mapping_generator import MappingGenerator
from os.path import realpath, dirname

file_dp = dirname(realpath(__file__))
__LIB_ROOT__DIR__ = dirname(file_dp)

__name__ = "biodblinker"
__description__ = "A library for linking entities of biological knowledge bases."
__version__ = "0.0.1"
__data_version__ = "0.0.1"
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
    "WikipediaLinker", "CellosaurusLinker", 'MappingGenerator'
    "__version__", "__name__", "__description__"
]

