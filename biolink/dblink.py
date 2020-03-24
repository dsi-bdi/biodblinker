from biolink.config import load_data_map, get_database_mapping_sources
from abc import ABC, abstractmethod


class DatabaseLinker(ABC):
    """ An abstract class for a database linker
    """
    def __init__(self):
        """
        Initialize a DatabaseLinker object
        """
        self._database_name = ""
        self._database_mapping_filenames = dict()
        self._database_linking_dictionaries = dict()
        self._database_linking_dictionaries_paths = dict()

    def _load_linking_dictionaries(self):
        """ Load database linking dictionaries
        """
        self._database_linking_dictionaries_paths = get_database_mapping_sources(self._database_name)

        for mapping_name, mapping_filepath in self._database_linking_dictionaries_paths.items():
            self._database_linking_dictionaries[mapping_name] = load_data_map(mapping_filepath)

    @property
    def mapping_file_paths(self):
        """ Get the list of mapping file paths associated with the database
        """
        return list(self._database_linking_dictionaries_paths.values())

    @property
    def mapping_dictionary_names(self):
        """ Get the list of the mapping dictionary names associated with the database
        """
        return list(self._database_linking_dictionaries_paths.keys())

    def _convert_ids_using_target_dictionary(self, ids_list, dictionary_name):
        """ Convert a list of ids using a corresponding target dictionary

        Parameters
        ----------
        ids_list : list
            a list of ids
        dictionary_name : str
            the name of the target dictionary

        Returns
        -------
        list
            a list of corresponding ids extracted from target dictionary
        """
        target_dict = self._database_linking_dictionaries[dictionary_name]
        outcome_ids = []
        for entry in ids_list:
            if entry in target_dict:
                outcome_ids.append(target_dict[entry])
            else:
                outcome_ids.append([])
        return outcome_ids

    def _get_keys_of_target_dictionary(self, dictionary_name):
        """ Get keys of  a database mapping dictionary

        Parameters
        ----------
        dictionary_name : str
            the dictionary name

        Returns
        -------
        list
            list of dictionary keys
        """
        return set(self._database_linking_dictionaries[dictionary_name].keys())


class KEGGLinker(DatabaseLinker):
    """ KEGG database linker class
    """

    def __init__(self):
        """ Initialize an object of the class KEGGLinker
        """
        super(KEGGLinker, self).__init__()
        self._database_name = "kegg"

        self._load_linking_dictionaries()

    def convert_geneid_to_uniprot(self, gene_list):
        """ Convert a list of KEGG gene ids to uniprot accessions

        Parameters
        ----------
        gene_list : list
            a list of KEGG gene ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of uniprot accession codes
        """
        return self._convert_ids_using_target_dictionary(gene_list, "gene_to_uniprot")

    def convert_geneid_to_names(self, gene_list):
        """ Convert a list of KEGG gene ids to gene names

        Parameters
        ----------
        gene_list : list
            a list of KEGG gene ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of gene names
        """
        return self._convert_ids_using_target_dictionary(gene_list, "gene_to_names")

    def convert_geneid_to_ensembl(self, gene_list):
        """ Convert a list of KEGG gene ids to ensembl gene ids

        Parameters
        ----------
        gene_list : list
            a list of KEGG gene ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of ensembl gene ids
        """
        return self._convert_ids_using_target_dictionary(gene_list, "gene_to_ensembl")

    def convert_drugid_to_chebi(self, drug_list):
        """ Convert a list of KEGG drug ids to chebi chemical ids

        Parameters
        ----------
        drug_list : list
            a list of KEGG drug ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of chebi chemical ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drug_to_chebi")

    def convert_drugid_to_chembl(self, drug_list):
        """ Convert a list of KEGG drug ids to chembl compound ids

        Parameters
        ----------
        drug_list : list
            a list of KEGG drug ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of chembl compound ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drug_to_chembl")

    def convert_drugid_to_hmdb(self, drug_list):
        """ Convert a list of KEGG drug ids to hmdb metabolite ids

        Parameters
        ----------
        drug_list : list
            a list of KEGG drug ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of hmdb metabolite ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drug_to_hmdb")

    def convert_drugid_to_hsdb(self, drug_list):
        """ Convert a list of KEGG drug ids to hsdb substance ids

        Parameters
        ----------
        drug_list : list
            a list of KEGG drug ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of hsdb substance ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drug_to_hsdb")

    def convert_drugid_to_knapsack(self, drug_list):
        """ Convert a list of KEGG drug ids to knapsack metabolite ids

        Parameters
        ----------
        drug_list : list
            a list of KEGG drug ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of knapsack metabolite ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drug_to_knapsack")

    def convert_drugid_to_ligandbox(self, drug_list):
        """ Convert a list of KEGG drug ids to ligandbox compound ids

        Parameters
        ----------
        drug_list : list
            a list of KEGG drug ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of ligandbox compound ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drug_to_ligandbox")

    def convert_drugid_to_massbank(self, drug_list):
        """ Convert a list of KEGG drug ids to massbank accessions ids

        Parameters
        ----------
        drug_list : list
            a list of KEGG drug ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of massbank accessions ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drug_to_massbank")

    def convert_drugid_to_names(self, drug_list):
        """ Convert a list of KEGG drug ids to names

        Parameters
        ----------
        drug_list : list
            a list of KEGG drug ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of names
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drug_to_names")

    def convert_drugid_to_nikkaji(self, drug_list):
        """ Convert a list of KEGG drug ids to nikkaji chemical ids

        Parameters
        ----------
        drug_list : list
            a list of KEGG drug ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of nikkaji chemical ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drug_to_nikkaji")

    def convert_drugid_to_pdb_ccd(self, drug_list):
        """ Convert a list of KEGG drug ids to pdb chemical component ids

        Parameters
        ----------
        drug_list : list
            a list of KEGG drug ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of chemical component ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drug_to_pdb-ccd")

    def convert_drugid_to_pubchem_substance(self, drug_list):
        """ Convert a list of KEGG drug ids to pubchem substance ids

        Parameters
        ----------
        drug_list : list
            a list of KEGG drug ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of pubchem substance ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drug_to_pubchem")

    def convert_drugid_to_sider(self, drug_list):
        """ Convert a list of KEGG drug ids to sider drug ids

        Parameters
        ----------
        drug_list : list
            a list of KEGG drug ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of sider drug ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drug_to_sider")

    def convert_compound_to_3dmet(self, compound_list):
        """ Convert a list of KEGG compound ids to 3dmet molecule ids

        Parameters
        ----------
        compound_list : list
            a list of KEGG compound ids e.g. ['C17612', 'C21587']

        Returns
        -------
        list
            a list of lists of 3dmet molecule ids
        """
        return self._convert_ids_using_target_dictionary(compound_list, "compound_to_3dmet")

    def convert_compound_to_chebi(self, compound_list):
        """ Convert a list of KEGG compound ids to chebi chemical ids

        Parameters
        ----------
        compound_list : list
            a list of KEGG compound ids e.g. ['C17612', 'C21587']

        Returns
        -------
        list
            a list of lists of chebi chemical ids
        """
        return self._convert_ids_using_target_dictionary(compound_list, "compound_to_chebi")

    def convert_compound_to_chembl(self, compound_list):
        """ Convert a list of KEGG compound ids to chembl compound ids

        Parameters
        ----------
        compound_list : list
            a list of KEGG compound ids e.g. ['C17612', 'C21587']

        Returns
        -------
        list
            a list of lists of chembl compound ids
        """
        return self._convert_ids_using_target_dictionary(compound_list, "compound_to_chembl")

    def convert_compound_to_hmdb(self, compound_list):
        """ Convert a list of KEGG compound ids to hmdb metabolite ids

        Parameters
        ----------
        compound_list : list
            a list of KEGG compound ids e.g. ['C17612', 'C21587']

        Returns
        -------
        list
            a list of lists of hmdb metabolite ids
        """
        return self._convert_ids_using_target_dictionary(compound_list, "compound_to_hmdb")

    def convert_compound_to_hsdb(self, compound_list):
        """ Convert a list of KEGG compound ids to hsdb substance ids

        Parameters
        ----------
        compound_list : list
            a list of KEGG compound ids e.g. ['C17612', 'C21587']

        Returns
        -------
        list
            a list of lists of hsdb substance ids
        """
        return self._convert_ids_using_target_dictionary(compound_list, "compound_to_hsdb")

    def convert_compound_to_knapsack(self, compound_list):
        """ Convert a list of KEGG compound ids to knapsack metabolite ids

        Parameters
        ----------
        compound_list : list
            a list of KEGG compound ids e.g. ['C17612', 'C21587']

        Returns
        -------
        list
            a list of lists of knapsack metabolite ids
        """
        return self._convert_ids_using_target_dictionary(compound_list, "compound_to_knapsack")

    def convert_compound_to_lipidbank(self, compound_list):
        """ Convert a list of KEGG compound ids to lipidbank lipid ids

        Parameters
        ----------
        compound_list : list
            a list of KEGG compound ids e.g. ['C17612', 'C21587']

        Returns
        -------
        list
            a list of lists of lipidbank lipid ids
        """
        return self._convert_ids_using_target_dictionary(compound_list, "compound_to_lipidbank")

    def convert_compound_to_lipidmaps(self, compound_list):
        """ Convert a list of KEGG compound ids to lipidmaps lipid ids

        Parameters
        ----------
        compound_list : list
            a list of KEGG compound ids e.g. ['C17612', 'C21587']

        Returns
        -------
        list
            a list of lists of lipidmaps lipid ids
        """
        return self._convert_ids_using_target_dictionary(compound_list, "compound_to_lipidmaps")

    def convert_compound_to_massbank(self, compound_list):
        """ Convert a list of KEGG compound ids to massbank accessions ids

        Parameters
        ----------
        compound_list : list
            a list of KEGG compound ids e.g. ['C17612', 'C21587']

        Returns
        -------
        list
            a list of lists of massbank accessions ids
        """
        return self._convert_ids_using_target_dictionary(compound_list, "compound_to_massbank")

    def convert_compound_to_names(self, compound_list):
        """ Convert a list of KEGG compound ids to names

        Parameters
        ----------
        compound_list : list
            a list of KEGG compound ids e.g. ['C17612', 'C21587']

        Returns
        -------
        list
            a list of lists of compound names
        """
        return self._convert_ids_using_target_dictionary(compound_list, "compound_to_names")

    def convert_compound_to_nikkaji(self, compound_list):
        """ Convert a list of KEGG compound ids to nikkaji chemical ids

        Parameters
        ----------
        compound_list : list
            a list of KEGG compound ids e.g. ['C17612', 'C21587']

        Returns
        -------
        list
            a list of lists of nikkaji chemical ids
        """
        return self._convert_ids_using_target_dictionary(compound_list, "compound_to_nikkaji")

    def convert_compound_to_pdb_ccd(self, compound_list):
        """ Convert a list of KEGG compound ids to pdb chemical component ids

        Parameters
        ----------
        compound_list : list
            a list of KEGG compound ids e.g. ['C17612', 'C21587']

        Returns
        -------
        list
            a list of lists of pdb chemical component ids
        """
        return self._convert_ids_using_target_dictionary(compound_list, "compound_to_pdb-ccd")

    def convert_compound_to_pubchem_substance(self, compound_list):
        """ Convert a list of KEGG compound ids to pubchem substance ids

        Parameters
        ----------
        compound_list : list
            a list of KEGG compound ids e.g. ['C17612', 'C21587']

        Returns
        -------
        list
            a list of lists of pubchem substance ids
        """
        return self._convert_ids_using_target_dictionary(compound_list, "compound_to_pubchem")

    def convert_brite_to_names(self, brite_list):
        """ Convert a list of KEGG brite ids to brite names

        Parameters
        ----------
        brite_list : list
            a list of KEGG brite ids e.g. ['br08901', 'br08902']

        Returns
        -------
        list
            a list of lists of brite names
        """
        return self._convert_ids_using_target_dictionary(brite_list, "brite_to_names")

    def convert_drug_group_to_names(self, drug_group_list):
        """ Convert a list of KEGG drug group ids to drug group names

        Parameters
        ----------
        drug_group_list : list
            a list of KEGG drug group ids e.g. ['DG00001', 'DG00002']

        Returns
        -------
        list
            a list of lists of drug group names
        """
        return self._convert_ids_using_target_dictionary(drug_group_list, "dgroup_to_names")

    def convert_disease_to_names(self, disease_list):
        """ Convert a list of KEGG disease ids to disease names

        Parameters
        ----------
        disease_list : list
            a list of KEGG disease ids e.g. ['H00001', 'H00002']

        Returns
        -------
        list
            a list of lists of disease names
        """
        return self._convert_ids_using_target_dictionary(disease_list, "disease_to_names")

    def convert_disease_to_omim(self, disease_list):
        """ Convert a list of KEGG disease ids to omim disease ids

        Parameters
        ----------
        disease_list : list
            a list of KEGG disease ids e.g. ['H00001', 'H00002']

        Returns
        -------
        list
            a list of lists of omim disease ids
        """
        return self._convert_ids_using_target_dictionary(disease_list, "disease_to_omim")

    def convert_disease_to_pubmed(self, disease_list):
        """ Convert a list of KEGG disease ids to pubmed article ids

        Parameters
        ----------
        disease_list : list
            a list of KEGG disease ids e.g. ['H00001', 'H00002']

        Returns
        -------
        list
            a list of lists of pubmed article ids
        """
        return self._convert_ids_using_target_dictionary(disease_list, "disease_to_pubmed")

    def convert_environ_to_names(self, environ_list):
        """ Convert a list of KEGG environ ids to environ names

        Parameters
        ----------
        environ_list : list
            a list of KEGG environ ids e.g. ['E00001', 'E00002']

        Returns
        -------
        list
            a list of lists of environ names
        """
        return self._convert_ids_using_target_dictionary(environ_list, "environ_to_names")

    def convert_enzyme_to_names(self, enzyme_list):
        """ Convert a list of KEGG enzyme ids to enzyme names

        Parameters
        ----------
        enzyme_list : list
            a list of KEGG enzyme ids e.g. ['ec:1.1.1.5', 'ec:1.1.1.7']

        Returns
        -------
        list
            a list of lists of enzyme names
        """
        return self._convert_ids_using_target_dictionary(enzyme_list, "enzyme_to_names")

    def convert_glycan_to_names(self, glycan_list):
        """ Convert a list of KEGG glycan ids to glycan names

        Parameters
        ----------
        glycan_list : list
            a list of KEGG glycan ids e.g. ['G00001', 'G00002']

        Returns
        -------
        list
            a list of lists of glycan namess
        """
        return self._convert_ids_using_target_dictionary(glycan_list, "glycan_to_names(")

    def convert_ko_to_names(self, ontology_list):
        """ Convert a list of KEGG ontology entry ids to ontology entry names

        Parameters
        ----------
        ontology_list : list
            a list of KEGG ontology entry ids e.g. ['K00001', 'K00002']

        Returns
        -------
        list
            a list of lists of ontology entry names
        """
        return self._convert_ids_using_target_dictionary(ontology_list, "ko_to_names")

    def convert_module_to_names(self, module_list):
        """ Convert a list of KEGG module ids to module names

        Parameters
        ----------
        module_list : list
            a list of KEGG module ids e.g. ['M00001', 'M00002']

        Returns
        -------
        list
            a list of lists of module namess
        """
        return self._convert_ids_using_target_dictionary(module_list, "module_to_names")

    def convert_network_to_names(self, network_list):
        """ Convert a list of KEGG network ids to network names

        Parameters
        ----------
        network_list : list
            a list of KEGG network ids e.g. ['N00001', 'N00002']

        Returns
        -------
        list
            a list of lists of network names
        """
        return self._convert_ids_using_target_dictionary(network_list, "network_to_names")

    def convert_pathway_to_names(self, pathway_list):
        """ Convert a list of KEGG pathway ids to pathway names

        Parameters
        ----------
        pathway_list : list
            a list of KEGG pathway ids e.g. ['map00010', 'map00020']

        Returns
        -------
        list
            a list of lists of pathway names
        """
        return self._convert_ids_using_target_dictionary(pathway_list, "pathway_to_names")

    def convert_reaction_to_names(self, reaction_list):
        """ Convert a list of KEGG reaction ids to reaction names

        Parameters
        ----------
        reaction_list : list
            a list of KEGG reaction ids e.g. ['R00001', 'R00002']

        Returns
        -------
        list
            a list of lists of reaction names
        """
        return self._convert_ids_using_target_dictionary(reaction_list, "reaction_to_names")

    def convert_variant_to_names(self, variant_list):
        """ Convert a list of KEGG variant ids to variant names

        Parameters
        ----------
        variant_list : list
            a list of KEGG variant ids e.g. ['hsa_var:100v1', 'hsa_var:1019v1']

        Returns
        -------
        list
            a list of lists of variant names
        """
        return self._convert_ids_using_target_dictionary(variant_list, "variant_to_names")
    @property
    def gene_ids(self):
        """
        All available gene ids in the database
        """
        return self._get_keys_of_target_dictionary("gene_to_names")

    @property
    def drug_ids(self):
        """
        All available drug ids in the database
        """
        return self._get_keys_of_target_dictionary("drug_to_names")

    @property
    def brite_ids(self):
        """
        All available brite ids in the database
        """
        return self._get_keys_of_target_dictionary("brite_to_names")

    @property
    def compound_ids(self):
        """
        All available compound ids in the database
        """
        return self._get_keys_of_target_dictionary("compound_to_names")

    @property
    def disease_ids(self):
        """
        All available disease ids in the database
        """
        return self._get_keys_of_target_dictionary("disease_to_names")

    @property
    def environ_ids(self):
        """
        All available environ ids in the database
        """
        return self._get_keys_of_target_dictionary("environ_to_names")

    @property
    def enzyme_ids(self):
        """
        All available enzyme ids in the database
        """
        return self._get_keys_of_target_dictionary("enzyme_to_names")

    @property
    def network_ids(self):
        """
        All available network ids in the database
        """
        return self._get_keys_of_target_dictionary("network_to_names")

    @property
    def pathway_ids(self):
        """
        All available pathway ids in the database
        """
        return self._get_keys_of_target_dictionary("pathway_to_names")

    @property
    def reaction_ids(self):
        """
        All available reaction ids in the database
        """
        return self._get_keys_of_target_dictionary("reaction_to_names")

    @property
    def variant_ids(self):
        """
        All available variant ids in the database
        """
        return self._get_keys_of_target_dictionary("variant_to_names")


class DrugBankLinker(DatabaseLinker):
    """ Drugbank database linker class
    """

    def __init__(self):
        """ Initialize an object of the class DrugBankLinker
        """
        super(DrugBankLinker, self).__init__()
        self._database_name = "drugbank"

        self._load_linking_dictionaries()

    def convert_drugs_to_bindingdb(self, drug_list):
        """ Convert a list of DrugBank drug ids to bindingdb compound ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of bindingdb compound ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_bindingdb")

    def convert_drugs_to_chebi(self, drug_list):
        """ Convert a list of DrugBank drug ids to chebi chemical ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of chebi chemical ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_chebi")

    def convert_drugs_to_chembl(self, drug_list):
        """ Convert a list of DrugBank drug ids to chembl compound ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of chembl compound ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_chembl")

    def convert_drugs_to_chemspider(self, drug_list):
        """ Convert a list of DrugBank drug ids to chemspider chemical ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of chemspider chemical ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_chemspider")

    def convert_drugs_to_dpd(self, drug_list):
        """ Convert a list of DrugBank drug ids to drug product database drug ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of drug product database drug ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_dpd")

    def convert_drugs_to_genbank(self, drug_list):
        """ Convert a list of DrugBank drug ids to genbank nucleotide ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of guide to genbank nucleotide ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_genbank")

    def convert_drugs_to_guide_pharmacology(self, drug_list):
        """ Convert a list of DrugBank drug ids to guide to pharmacology ligand ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of guide to pharmacology ligand ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_pharma_guide")

    def convert_drugs_to_iuphar(self, drug_list):
        """ Convert a list of DrugBank drug ids to iuphar ligand ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of iuphar ligand ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_iuphar")

    def convert_drugs_to_kegg_compound(self, drug_list):
        """ Convert a list of DrugBank drug ids to kegg compound ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of kegg compound ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_kegg_compound")

    def convert_drugs_to_kegg_drug(self, drug_list):
        """ Convert a list of DrugBank drug ids to kegg drug ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of kegg drug ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_kegg_drug")

    def convert_drugs_to_pdb(self, drug_list):
        """ Convert a list of DrugBank drug ids to protein data bank ligand ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of protein data bank ligand ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_pdb")

    def convert_drugs_to_pharmgkb(self, drug_list):
        """ Convert a list of DrugBank drug ids to pharmgkb mocule ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of pharmgkb mocule ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_pharmgkb")

    def convert_drugs_to_pubchem_compound(self, drug_list):
        """ Convert a list of DrugBank drug ids to pubchem compound ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of pubchem compound ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_pubchem_compound")

    def convert_drugs_to_pubchem_substance(self, drug_list):
        """ Convert a list of DrugBank drug ids to pubchem substance ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of pubchem substance ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_pubchem_substance")

    def convert_drugs_to_ttd(self, drug_list):
        """ Convert a list of DrugBank drug ids to therapeutic targets database drug ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of therapeutic targets database drug ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_ttd")

    def convert_drugs_to_uniprotkb(self, drug_list):
        """ Convert a list of DrugBank drug ids to uniprotkb ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of uniprotkb ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_uniprotkb")

    def convert_drugs_to_wikipedia(self, drug_list):
        """ Convert a list of DrugBank drug ids to wikipedia entries

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of wikipedia entries
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_wikipedia")

    def convert_drugs_to_secondaryids(self, drug_list):
        """ Convert a list of DrugBank drug ids to secondary drug ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of secondary drug ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_ids")

    def convert_drugs_to_names(self, drug_list):
        """ Convert a list of DrugBank drug ids to drug names

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of drug names
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_names")

    @property
    def drug_ids(self):
        """
        All available drug ids in the database
        """
        return self._get_keys_of_target_dictionary("drugbank_to_names")


class SiderLinker(DatabaseLinker):
    """ Sider database linker class
    """

    def __init__(self):
        """ Initialize an object of the class SiderLinker
        """
        super(SiderLinker, self).__init__()
        self._database_name = "sider"

        self._load_linking_dictionaries()

    def convert_drugs_to_bindingdb(self, drug_list):
        """ Convert a list of Sider drug ids to kegg drug ids

        Parameters
        ----------
        drug_list : list
            a list of Sider drug ids e.g. ['CID100000961', 'CID100000206']

        Returns
        -------
        list
            a list of lists of kegg drug ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "sider_to_bindingdb")

    def convert_drugs_to_chebi(self, drug_list):
        """ Convert a list of Sider drug ids to chebi chemical ids

        Parameters
        ----------
        drug_list : list
            a list of Sider drug ids e.g. ['CID100000961', 'CID100000206']

        Returns
        -------
        list
            a list of lists of chebi chemical ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "sider_to_chebi")

    def convert_drugs_to_chembl(self, drug_list):
        """ Convert a list of Sider drug ids to chembl compound ids

        Parameters
        ----------
        drug_list : list
            a list of Sider drug ids e.g. ['CID100000961', 'CID100000206']

        Returns
        -------
        list
            a list of lists of chembl compound ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "sider_to_chembl")

    def convert_drugs_to_drugbank(self, drug_list):
        """ Convert a list of Sider drug ids to drugbank drug ids

        Parameters
        ----------
        drug_list : list
            a list of Sider drug ids e.g. ['CID100000961', 'CID100000206']

        Returns
        -------
        list
            a list of lists of drugbank drug ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "sider_to_drugbank")

    def convert_drugs_to_kegg_drug(self, drug_list):
        """ Convert a list of Sider drug ids to kegg drug ids

        Parameters
        ----------
        drug_list : list
            a list of Sider drug ids e.g. ['CID100000961', 'CID100000206']

        Returns
        -------
        list
            a list of lists of kegg drug ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "sider_to_kegg")

    def convert_drugs_to_pubchem_substance(self, drug_list):
        """ Convert a list of Sider drug ids to pubchem substance ids

        Parameters
        ----------
        drug_list : list
            a list of Sider drug ids e.g. ['CID100000961', 'CID100000206']

        Returns
        -------
        list
            a list of lists of pubchem substance ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "sider_to_pubchem_substance")

    def convert_drugs_to_pubchem_compound(self, drug_list):
        """ Convert a list of Sider drug ids to pubchem compound ids

        Parameters
        ----------
        drug_list : list
            a list of Sider drug ids e.g. ['CID100000961', 'CID100000206']

        Returns
        -------
        list
            a list of lists of pubchem compound ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "sider_to_pubchem_compound")

    def convert_drugs_to_atc(self, drug_list):
        """ Convert a list of Sider drug ids to atc codes

        Parameters
        ----------
        drug_list : list
            a list of Sider drug ids e.g. ['CID100000961', 'CID100000206']

        Returns
        -------
        list
            a list of lists of atc codes
        """
        return self._convert_ids_using_target_dictionary(drug_list, "sider_to_atc")

    def convert_drugs_to_names(self, drug_list):
        """ Convert a list of Sider drug ids to drug names

        Parameters
        ----------
        drug_list : list
            a list of Sider drug ids e.g. ['CID100000961', 'CID100000206']

        Returns
        -------
        list
            a list of lists of drug names
        """
        return self._convert_ids_using_target_dictionary(drug_list, "sider_to_names")

    @property
    def drug_ids(self):
        """
        All available drug ids in the database
        """
        return self._get_keys_of_target_dictionary("sider_to_names")


class UniprotLinker(DatabaseLinker):
    """ Uniprot database linker class
    """

    def __init__(self):
        """ Initialize an object of the class UniprotLinker
        """
        super(UniprotLinker, self).__init__()
        self._database_name = "uniprot"

        self._load_linking_dictionaries()

    def convert_uniprot_to_biogrid(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to biogrid protein ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of biogrid protein ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_biogrid")

    def convert_uniprot_to_chembl(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to chembl compound ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of chembl compound ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_chembl")

    def convert_uniprot_to_dip(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to dip protein ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of dip protein ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_dip")

    def convert_uniprot_to_embl_cds(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to embl-cds sequences

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of embl-cds sequences
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_embl-cds")

    def convert_uniprot_to_embl(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to embl sequences

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of embl sequences
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_embl")

    def convert_uniprot_to_ensembl(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to ensembl gene ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of ensembl gene ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_ensembl")

    def convert_uniprot_to_genedb(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to genedb gene ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of genedb gene ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_genedb")

    def convert_uniprot_to_gene_name(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to gene names

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of biogrid gene names
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_gene_name")

    def convert_uniprot_to_gene_synonym(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to gene synonyms

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of gene synonyms
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_gene_synonym")

    def convert_uniprot_to_gi(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to genebank gene ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of genebank gene ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_gi")

    def convert_uniprot_to_hgnc(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to hgnc gene ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of hgnc gene ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_hgnc")

    def convert_uniprot_to_hpa(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to hpa protein ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of hpa protein ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_hpa")

    def convert_uniprot_to_ko(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to kegg ontology ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of kegg ontology ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_ko")

    def convert_uniprot_to_kegg(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to kegg gene ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of kegg gene ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_kegg")

    def convert_uniprot_to_mim(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to omim entry ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of omim entry ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_mim")

    def convert_uniprot_to_mint(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to mint protein ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of mint protein ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_mint")

    def convert_uniprot_to_pharmgkb(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to pharmkgb gene ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of pharmkgb gene ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_pharmgkb")

    def convert_uniprot_to_proteomicsdb(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to proteomicsdb protein ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of proteomicsdb protein ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_proteomicsdb")

    def convert_uniprot_to_refseq(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to refseq protein ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of refseq protein ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_refseq")

    def convert_uniprot_to_string(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to string protein ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of string protein ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_string")

    def convert_uniprot_to_uniparc(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to uniparc protein ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of uniparc protein ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_uniparc")

    def convert_uniprot_to_uniprotkb_id(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to uniprotkb protein ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of uniprotkb protein ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_uniprotkb-id")

    def convert_uniprot_to_uniref100(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to uniref100 protein ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of uniref100 protein ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_uniref100")

    def convert_uniprot_to_names(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to protein names

        Parameters
        ----------
        acc_list : list
            a list of Uniprot drug ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of protein names
        """
        return self._convert_ids_using_target_dictionary(acc_list, "uniprot_names")

    @property
    def uniprot_ids(self):
        """
        All available accession ids in the database
        """
        return self._get_keys_of_target_dictionary("uniprot_names")

class GeneNameLinker(DatabaseLinker):
    """ Gene name database linker class
    """

    def __init__(self):
        """ Initialize an object of the class GeneNameLinker
        """
        super(GeneNameLinker, self).__init__()
        self._database_name = "gene_name"

        self._load_linking_dictionaries()

    def convert_gene_name_to_uniprot(self, gene_name_list):
        """ Convert a list of gene names to uniprot accessions

        Parameters
        ----------
        gene_name_list : list
            a list of gene names e.g. ['LATS1', 'COX2']

        Returns
        -------
        list
            a list of uniprot accession codes
        """
        return self._convert_ids_using_target_dictionary(gene_name_list, "gene_name_to_uniprot")
