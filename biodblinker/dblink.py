from biodblinker.config import load_data_map, get_database_mapping_sources, verify_biodblinker
from abc import ABC, abstractmethod


class DatabaseLinker(ABC):
    """ An abstract class for a database linker
    """
    def __init__(self):
        """
        Initialize a DatabaseLinker object
        """
        verify_biodblinker()
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

    def convert_names_to_geneids(self, name_list):
        """ Convert a list of KEGG gene namse to gene ids

        Parameters
        ----------
        name_list : list
            a list of KEGG gene ids e.g. ['GRF2', 'ND6']

        Returns
        -------
        list
            a list of lists of kegg gene
        """
        return self._convert_ids_using_target_dictionary(name_list, "names_to_genes")

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

    def convert_names_to_drugids(self, name_list):
        """ Convert a list of KEGG drug names to drug ids

        Parameters
        ----------
        name_list : list
            a list of KEGG drug names e.g. ['Danicamtiv', 'Balstilimab']

        Returns
        -------
        list
            a list of lists of drug ids
        """
        return self._convert_ids_using_target_dictionary(name_list, "names_to_drugs")

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

    def convert_drugid_to_drugbank(self, drug_list):
        """ Convert a list of KEGG drug ids to drugbank drug ids

        Parameters
        ----------
        drug_list : list
            a list of KEGG drug ids e.g. ['D07630', 'D02238']

        Returns
        -------
        list
            a list of lists of drugbank drug ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drug_to_drugbank")

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

    def convert_names_to_diseases(self, name_list):
        """ Convert a list of KEGG disease names to disease ids

        Parameters
        ----------
        name_list : list
            a list of KEGG disease names e.g. ['Trichosporonosis', 'Erdheim-Chester disease']

        Returns
        -------
        list
            a list of lists of disease ids
        """
        return self._convert_ids_using_target_dictionary(name_list, "names_to_disease")

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

    def convert_disease_to_mesh(self, disease_list):
        """ Convert a list of KEGG disease ids to mesh disease ids

        Parameters
        ----------
        disease_list : list
            a list of KEGG disease ids e.g. ['H00001', 'H00002']

        Returns
        -------
        list
            a list of lists of mesh disease ids
        """
        return self._convert_ids_using_target_dictionary(disease_list, "disease_to_mesh")

    def convert_disease_to_icd_10(self, disease_list):
        """ Convert a list of KEGG disease ids to icd 10 disease ids

        Parameters
        ----------
        disease_list : list
            a list of KEGG disease ids e.g. ['H00001', 'H00002']

        Returns
        -------
        list
            a list of lists of icd 10 disease ids
        """
        return self._convert_ids_using_target_dictionary(disease_list, "disease_to_icd_10")
    
    def convert_disease_to_icd_11(self, disease_list):
        """ Convert a list of KEGG disease ids to icd 11 disease ids

        Parameters
        ----------
        disease_list : list
            a list of KEGG disease ids e.g. ['H00001', 'H00002']

        Returns
        -------
        list
            a list of lists of icd 11 disease ids
        """
        return self._convert_ids_using_target_dictionary(disease_list, "disease_to_icd_11")

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

    def convert_glycan_to_names(self, glycan_list):
        """ Convert a list of KEGG glycan ids to glycan names

        Parameters
        ----------
        glycan_list : list
            a list of KEGG glycan ids e.g. ['G00001', 'G00002']

        Returns
        -------
        list
            a list of lists of glycan names
        """
        return self._convert_ids_using_target_dictionary(glycan_list, "glycan_to_names")

    def convert_names_to_glycans(self, name_list):
        """ Convert a list of KEGG glycan names to glycan ids

        Parameters
        ----------
        name_list : list
            a list of KEGG glycan names e.g. ['(GlcNAc)2 (Man)1 (PP-Dol)1', 'N-Acetyl-D-glucosaminyldiphosphodolichol']

        Returns
        -------
        list
            a list of lists of glycan ids
        """
        return self._convert_ids_using_target_dictionary(name_list, "names_to_glycans")

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

    def convert_names_to_networks(self, name_list):
        """ Convert a list of KEGG network names to network ids

        Parameters
        ----------
        name_list : list
            a list of KEGG network names e.g. ['Hepatitis B virus (HBV)', 'Epstein-Barr virus (EBV)']

        Returns
        -------
        list
            a list of lists of network ids
        """
        return self._convert_ids_using_target_dictionary(name_list, "names_to_networks")

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

    def convert_names_to_pathways(self, name_list):
        """ Convert a list of KEGG pathway names to pathway ids

        Parameters
        ----------
        name_list : list
            a list of KEGG pathway ids e.g. ['Galactose metabolism', 'Fatty acid biosynthesis']

        Returns
        -------
        list
            a list of lists of pathway ids
        """
        return self._convert_ids_using_target_dictionary(name_list, "names_to_pathways")

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
    def disease_ids(self):
        """
        All available disease ids in the database
        """
        return self._get_keys_of_target_dictionary("disease_to_names")

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
        """ Convert a list of DrugBank drug ids to pharmgkb molecule ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of pharmgkb molecule ids
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

    def convert_drugs_to_uniprot(self, drug_list):
        """ Convert a list of DrugBank drug ids to uniprot accession ids

        Parameters
        ----------
        drug_list : list
            a list of DrugBank drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of uniprot accession ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "drugbank_to_uniprot")

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

    def convert_secondaryids_to_drugs(self, secondary_list):
        """ Convert a list of DrugBank secondary drug ids to drug ids

        Parameters
        ----------
        secondary_list : list
            a list of DrugBank secondary_list drug ids e.g. ['DB00001', 'DB00002']

        Returns
        -------
        list
            a list of lists of drug ids
        """
        return self._convert_ids_using_target_dictionary(secondary_list, "ids_to_drugbank")

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

    def convert_names_to_drugs(self, name_list):
        """ Convert a list of DrugBank drug names to drug ids

        Parameters
        ----------
        name_list : list
            a list of DrugBank drug ids e.g. ['Lepirudin', 'Cetuximab']

        Returns
        -------
        list
            a list of lists of drug ids
        """
        return self._convert_ids_using_target_dictionary(name_list, "names_to_drugbank")

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
        """ Convert a list of Sider drug ids to bindingdb drug ids

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

    def convert_names_to_drugs(self, name_list):
        """ Convert a list of Sider drug names to drug ids

        Parameters
        ----------
        name_list : list
            a list of Sider drug names e.g. ['carnitine', 'prostacyclin']

        Returns
        -------
        list
            a list of lists of drug ids
        """
        return self._convert_ids_using_target_dictionary(drug_list, "names_to_sider")

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of dip protein ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_dip")

    def convert_uniprot_to_drugbank(self, protein_list):
        """ Convert a list of Uniprot protein ids to Brugbank drug ids

        Parameters
        ----------
        protein_list : list
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of Drugbank drug ids
        """
        return self._convert_ids_using_target_dictionary(protein_list, "acc_to_drugbank")

    def convert_uniprot_to_embl_cds(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to embl-cds sequences

        Parameters
        ----------
        acc_list : list
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of hpa protein ids
        """
        return self._convert_ids_using_target_dictionary(acc_list, "acc_to_hpa")

    def convert_uniprot_to_kegg(self, acc_list):
        """ Convert a list of Uniprot protein accession ids to kegg gene ids

        Parameters
        ----------
        acc_list : list
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

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
            a list of Uniprot protein ids e.g. ['Q7L9L4', 'P42768']

        Returns
        -------
        list
            a list of lists of protein names
        """
        return self._convert_ids_using_target_dictionary(acc_list, "uniprot_names")

    def convert_names_to_uniprot(self, name_list):
        """ Convert a list of protein names to Uniprot protein accession ids

        Parameters
        ----------
        name_list : list
            a list of Uniprot protein ids e.g. ['001R', '003L']

        Returns
        -------
        list
            a list of lists of protein accessions
        """
        return self._convert_ids_using_target_dictionary(name_list, "names_to_uniprot")

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


class ChemblLinker(DatabaseLinker):
    """ Chembl database linker class
    """

    def __init__(self):
        """ Initialize an object of the class ChemblLinker
        """
        super(ChemblLinker, self).__init__()
        self._database_name = "chembl"

        self._load_linking_dictionaries()


    def convert_compound_to_drugbank_drugs(self, compound_list):
        """ Convert a list of Chembl compound ids to drugbank drugs

        Parameters
        ----------
        compound_list : list
            a list of chembl compound ids e.g. ['CHEMBL3710408', 'CHEMBL1293296']

        Returns
        -------
        list
            a list of lists of drugbank drugs
        """
        return self._convert_ids_using_target_dictionary(compound_list, "chembl_to_drugbank")

    def convert_compound_to_kegg_drugs(self, compound_list):
        """ Convert a list of Chembl compound ids to kegg drugs

        Parameters
        ----------
        compound_list : list
            a list of chembl compound ids e.g. ['CHEMBL3710408', 'CHEMBL1293296']

        Returns
        -------
        list
            a list of lists of kegg drugs
        """
        return self._convert_ids_using_target_dictionary(compound_list, "chembl_to_kegg_drug")
    
    def convert_compound_to_uniprot(self, compound_list):
        """ Convert a list of Chembl compound ids to uniprot accessions

        Parameters
        ----------
        compound_list : list
            a list of chembl compound ids e.g. ['CHEMBL3710408', 'CHEMBL1293296']

        Returns
        -------
        list
            a list of lists of uniprot accession codes
        """
        return self._convert_ids_using_target_dictionary(compound_list, "chembl_to_uniprot")



class BiogridLinker(DatabaseLinker):
    """ Biogrid database linker class
    """

    def __init__(self):
        """ Initialize an object of the class BiogridLinker
        """
        super(BiogridLinker, self).__init__()
        self._database_name = "biogrid"

        self._load_linking_dictionaries()

    def convert_biogrid_to_uniprot(self, protein_list):
        """ Convert a list of biogrid protein ids to Uniprot protein accession ids

        Parameters
        ----------
        protein_list : list
            a list of biogrid protein ids e.g. ['113361', '113365']

        Returns
        -------
        list
            a list of lists of uniprot protein accessions
        """
        return self._convert_ids_using_target_dictionary(protein_list, "biogrid_to_uniprot")


class DipLinker(DatabaseLinker):
    """ Dip database linker class
    """

    def __init__(self):
        """ Initialize an object of the class DipLinker
        """
        super(DipLinker, self).__init__()
        self._database_name = "dip"

        self._load_linking_dictionaries()

    def convert_dip_to_uniprot(self, protein_list):
        """ Convert a list of dip protein ids to Uniprot protein accession ids

        Parameters
        ----------
        protein_list : list
            a list of dip protein ids e.g. ['DIP-36676N', 'DIP-27584N']

        Returns
        -------
        list
            a list of lists of uniprot protein accessions
        """
        return self._convert_ids_using_target_dictionary(protein_list, "dip_to_uniprot")


class EmblLinker(DatabaseLinker):
    """ Embl database linker class
    """

    def __init__(self):
        """ Initialize an object of the class EmblLinker
        """
        super(EmblLinker, self).__init__()
        self._database_name = "embl"

        self._load_linking_dictionaries()

    def convert_embl_to_uniprot(self, sequence_list):
        """ Convert a list of embl sequence ids to Uniprot protein accession ids

        Parameters
        ----------
        sequence_list : list
            a list of embl sequence ids e.g. ['AK292717', 'AL008725']

        Returns
        -------
        list
            a list of lists of uniprot protein accessions
        """
        return self._convert_ids_using_target_dictionary(sequence_list, "embl_to_uniprot")


class EmblCDSLinker(DatabaseLinker):
    """ Embl-CDS database linker class
    """

    def __init__(self):
        """ Initialize an object of the class EmblCDSLinker
        """
        super(EmblCDSLinker, self).__init__()
        self._database_name = "embl-cds"

        self._load_linking_dictionaries()

    def convert_embl_to_uniprot(self, sequence_list):
        """ Convert a list of embl-cds sequence ids to Uniprot protein accession ids

        Parameters
        ----------
        sequence_list : list
            a list of embl-cds sequence ids e.g. ['EAW75893.1', 'AAC50710.1']

        Returns
        -------
        list
            a list of lists of uniprot protein accessions
        """
        return self._convert_ids_using_target_dictionary(sequence_list, "embl-cds_to_uniprot")


class EnsemblLinker(DatabaseLinker):
    """ Ensembl database linker class
    """

    def __init__(self):
        """ Initialize an object of the class EnsemblLinker
        """
        super(EnsemblLinker, self).__init__()
        self._database_name = "ensembl"

        self._load_linking_dictionaries()

    def convert_ensembl_to_uniprot(self, gene_list):
        """ Convert a list of ensembl ids to Uniprot protein accession ids

        Parameters
        ----------
        gene_list : list
            a list of ensembl gene ids e.g. ['EAW75893.1', 'AAC50710.1']

        Returns
        -------
        list
            a list of lists of uniprot protein accessions
        """
        return self._convert_ids_using_target_dictionary(gene_list, "ensembl_to_uniprot")

    def convert_ensembl_to_kegg_gene(self, gene_list):
        """ Convert a list of ensembl ids to Uniprot protein accession ids

        Parameters
        ----------
        gene_list : list
            a list of ensembl gene ids e.g. ['ENSG00000166913', 'ENSG00000175793']

        Returns
        -------
        list
            a list of lists of kegg gene ids
        """
        return self._convert_ids_using_target_dictionary(gene_list, "ensembl_to_kegg_gene")


class GeneDBLinker(DatabaseLinker):
    """ GeneDB database linker class
    """

    def __init__(self):
        """ Initialize an object of the class GeneDBLinker
        """
        super(GeneDBLinker, self).__init__()
        self._database_name = "genedb"

        self._load_linking_dictionaries()

    def convert_ensembl_to_uniprot(self, gene_list):
        """ Convert a list of genedb gene ids to Uniprot protein accession ids

        Parameters
        ----------
        gene_list : list
            a list of genedb gene ids e.g. ['Smp_196150.1:pep', 'Smp_009580.1:pep']

        Returns
        -------
        list
            a list of lists of uniprot protein accessions
        """
        return self._convert_ids_using_target_dictionary(gene_list, "genedb_to_uniprot")


class GILinker(DatabaseLinker):
    """ GI database linker class
    """

    def __init__(self):
        """ Initialize an object of the class GILinker
        """
        super(GILinker, self).__init__()
        self._database_name = "gi"

        self._load_linking_dictionaries()

    def convert_gi_to_uniprot(self, gene_list):
        """ Convert a list of gi gene ids to Uniprot protein accession ids

        Parameters
        ----------
        gene_list : list
            a list of gi gene ids e.g. ['377656702', '67464627']

        Returns
        -------
        list
            a list of lists of uniprot protein accessions
        """
        return self._convert_ids_using_target_dictionary(gene_list, "gi_to_uniprot")


class HGNCLinker(DatabaseLinker):
    """ HGNC database linker class
    """

    def __init__(self):
        """ Initialize an object of the class HGNCLinker
        """
        super(HGNCLinker, self).__init__()
        self._database_name = "hgnc"

        self._load_linking_dictionaries()

    def convert_hgnc_to_uniprot(self, gene_list):
        """ Convert a list of hgnc gene ids to Uniprot protein accession ids

        Parameters
        ----------
        gene_list : list
            a list of hgnc gene ids e.g. ['HGNC:12849', 'HGNC:12852']

        Returns
        -------
        list
            a list of lists of uniprot protein accessions
        """
        return self._convert_ids_using_target_dictionary(gene_list, "hgnc_to_uniprot")


class HPALinker(DatabaseLinker):
    """ HPA database linker class
    """

    def __init__(self):
        """ Initialize an object of the class HPALinker
        """
        super(HPALinker, self).__init__()
        self._database_name = "hpa"

        self._load_linking_dictionaries()

    def convert_hpa_to_uniprot(self, protein_list):
        """ Convert a list of hpa protein ids to Uniprot protein accession ids

        Parameters
        ----------
        protein_list : list
            a list of gi gene ids e.g. ['HPA011212', 'CAB016200']

        Returns
        -------
        list
            a list of lists of uniprot protein accessions
        """
        return self._convert_ids_using_target_dictionary(protein_list, "hpa_to_uniprot")


class MimLinker(DatabaseLinker):
    """ Mim database linker class
    """

    def __init__(self):
        """ Initialize an object of the class MimLinker
        """
        super(MimLinker, self).__init__()
        self._database_name = "mim"

        self._load_linking_dictionaries()

    def convert_mim_to_uniprot(self, entry_list):
        """ Convert a list of mim entry ids to Uniprot protein accession ids

        Parameters
        ----------
        entry_list : list
            a list of mim entry ids e.g. ['605066', '601288']

        Returns
        -------
        list
            a list of lists of uniprot protein accessions
        """
        return self._convert_ids_using_target_dictionary(entry_list, "mim_to_uniprot")

    def convert_mim_to_kegg_disease(self, entry_list):
        """ Convert a list of mim entry ids to KEGG disease ids

        Parameters
        ----------
        entry_list : list
            a list of mim entry ids e.g. ['605066', '601288']

        Returns
        -------
        list
            a list of lists of KEGG disease ids
        """
        return self._convert_ids_using_target_dictionary(entry_list, "mim_to_kegg_disease")

    def convert_mim_to_mesh(self, entry_list):
        """ Convert a list of mim entry ids to MeSH disease ids"

        Parameters
        ----------
        entry_list : list
            a list of mim entry ids e.g. ['605066', '601288']

        Returns
        -------
        list
            a list of lists of mesh disease ids
        """
        return self._convert_ids_using_target_dictionary(entry_list, "mim_to_mesh")


class MintLinker(DatabaseLinker):
    """ Mint database linker class
    """

    def __init__(self):
        """ Initialize an object of the class MintLinker
        """
        super(MintLinker, self).__init__()
        self._database_name = "mint"

        self._load_linking_dictionaries()

    def convert_mint_to_uniprot(self, protein_list):
        """ Convert a list of mint protein ids to Uniprot protein accession ids

        Parameters
        ----------
        protein_list : list
            a list of mint protein ids e.g. ['P62258', 'P31946']

        Returns
        -------
        list
            a list of lists of uniprot protein accessions
        """
        return self._convert_ids_using_target_dictionary(protein_list, "mint_to_uniprot")


class PharmgkbLinker(DatabaseLinker):
    """ PharmgkbLinker database linker class
    """

    def __init__(self):
        """ Initialize an object of the class PharmgkbLinker
        """
        super(PharmgkbLinker, self).__init__()
        self._database_name = "pharmgkb"

        self._load_linking_dictionaries()

    def convert_pharmgkb_to_uniprot(self, gene_list):
        """ Convert a list of pharmkgb gene ids to Uniprot protein accession ids

        Parameters
        ----------
        gene_list : list
            a list of pharmkgb gene ids e.g. ['PA37440', 'PA37444']

        Returns
        -------
        list
            a list of lists of uniprot protein accessions
        """
        return self._convert_ids_using_target_dictionary(gene_list, "pharmgkb_to_uniprot")

    def convert_pharmgkb_to_drugbank_drugs(self, molecule_list):
        """ Convert a list of pharmkgb molecule ids to Drugbank drug ids

        Parameters
        ----------
        molecule_list : list
            a list of pharmkgb molecule ids e.g. ['PA10318', 'PA10032']

        Returns
        -------
        list
            a list of lists of uniprot protein accessions
        """
        return self._convert_ids_using_target_dictionary(molecule_list, "pharmgkb_to_drugbank")


class ProteomicsDBLinker(DatabaseLinker):
    """ProteomicsDB database linker class
    """

    def __init__(self):
        """ Initialize an object of the class ProteomicsDBLinker
        """
        super(ProteomicsDBLinker, self).__init__()
        self._database_name = "proteomicsdb"

        self._load_linking_dictionaries()

    def convert_proteomicsdb_to_uniprot(self, protein_list):
        """ Convert a list of proteomicsdb protein ids to Uniprot protein accession ids

        Parameters
        ----------
        protein_list : list
            a list of proteomicsdb protein ids e.g. ['57377', '57378']

        Returns
        -------
        list
            a list of lists of Uniprot protein accession ids
        """
        return self._convert_ids_using_target_dictionary(protein_list, "proteomicsdb_to_uniprot")


class RefSeqLinker(DatabaseLinker):
    """RefSeq database linker class
    """

    def __init__(self):
        """ Initialize an object of the class RefSeqLinker
        """
        super(RefSeqLinker, self).__init__()
        self._database_name = "refseq"

        self._load_linking_dictionaries()

    def convert_refseq_to_uniprot(self, protein_list):
        """ Convert a list of refseq protein ids to Uniprot protein accession ids

        Parameters
        ----------
        protein_list : list
            a list of refseq protein ids e.g. ['NP_003395.1', 'NP_036611.2']

        Returns
        -------
        list
            a list of lists of Uniprot protein accession ids
        """
        return self._convert_ids_using_target_dictionary(protein_list, "refseq_to_uniprot")


class StringLinker(DatabaseLinker):
    """String database linker class
    """

    def __init__(self):
        """ Initialize an object of the class StringLinker
        """
        super(StringLinker, self).__init__()
        self._database_name = "string"

        self._load_linking_dictionaries()

    def convert_string_to_uniprot(self, protein_list):
        """ Convert a list of string protein ids to Uniprot protein accession ids

        Parameters
        ----------
        protein_list : list
            a list of string protein ids e.g. ['9606.ENSP00000361930', '9606.ENSP00000371267']

        Returns
        -------
        list
            a list of lists of Uniprot protein accession ids
        """
        return self._convert_ids_using_target_dictionary(protein_list, "string_to_uniprot")


class UniparcLinker(DatabaseLinker):
    """Uniparc database linker class
    """

    def __init__(self):
        """ Initialize an object of the class UniparcLinker
        """
        super(UniparcLinker, self).__init__()
        self._database_name = "uniparc"

        self._load_linking_dictionaries()

    def convert_uniparc_to_uniprot(self, protein_list):
        """ Convert a list of uniparc protein ids to Uniprot protein accession ids

        Parameters
        ----------
        protein_list : list
            a list of uniparc protein ids e.g. ['UPI000059C8F6', 'UPI000013CC64']

        Returns
        -------
        list
            a list of lists of Uniprot protein accession ids
        """
        return self._convert_ids_using_target_dictionary(protein_list, "uniparc_to_uniprot")


class UniprotKBLinker(DatabaseLinker):
    """UniprotKB database linker class
    """

    def __init__(self):
        """ Initialize an object of the class UniprotKBLinker
        """
        super(UniprotKBLinker, self).__init__()
        self._database_name = "uniprotkb"

        self._load_linking_dictionaries()

    def convert_uniprotkb_to_uniprot(self, protein_list):
        """ Convert a list of UniprotKB protein ids to Uniprot protein accession ids

        Parameters
        ----------
        protein_list : list
            a list of UniprotKB protein ids e.g. ['1433E_HUMAN', '1A1L1_HUMAN']

        Returns
        -------
        list
            a list of lists of Uniprot protein accession ids
        """
        return self._convert_ids_using_target_dictionary(protein_list, "uniprotkb_to_uniprot")


class Uniref100Linker(DatabaseLinker):
    """UniRef100 database linker class
    """

    def __init__(self):
        """ Initialize an object of the class Uniref100Linker
        """
        super(Uniref100Linker, self).__init__()
        self._database_name = "uniref100"

        self._load_linking_dictionaries()

    def convert_uniref100_to_uniprot(self, protein_list):
        """ Convert a list of UniRef100 protein ids to Uniprot protein accession ids

        Parameters
        ----------
        protein_list : list
            a list of UniRef100 protein ids e.g. ['9606.ENSP00000361930', '9606.ENSP00000371267']

        Returns
        -------
        list
            a list of lists of Uniprot protein accession ids
        """
        return self._convert_ids_using_target_dictionary(protein_list, "uniref100_to_uniprot")


class ChemspiderLinker(DatabaseLinker):
    """Chemspider database linker class
    """

    def __init__(self):
        """ Initialize an object of the class ChemspiderLinker
        """
        super(ChemspiderLinker, self).__init__()
        self._database_name = "chemspider"

        self._load_linking_dictionaries()

    def convert_chemspider_to_drugbank(self, molecule_list):
        """ Convert a list of chemspider chemical ids to DrugBank drug ids

        Parameters
        ----------
        molecule_list : list
            a list of chemspider chemical ids e.g. ['10482069', '4470656']

        Returns
        -------
        list
            a list of lists of DrugBank drug ids
        """
        return self._convert_ids_using_target_dictionary(molecule_list, "chemspider_to_drugbank")


class HMDBLinker(DatabaseLinker):
    """HMDB database linker class
    """

    def __init__(self):
        """ Initialize an object of the class HMDBLinker
        """
        super(HMDBLinker, self).__init__()
        self._database_name = "hmdb"

        self._load_linking_dictionaries()

    def convert_hmdb_to_kegg_drug(self, metabolite_list):
        """ Convert a list of hmdb metabolite ids to KEGG drug ids

        Parameters
        ----------
        metabolite_list : list
            a list of hmdb metabolite ids e.g. ['HMDB05032', 'HMDB41011']

        Returns
        -------
        list
            a list of lists of KEGG drug ids
        """
        return self._convert_ids_using_target_dictionary(metabolite_list, "hmdb_to_kegg_drug")


class GuidePharmaLinker(DatabaseLinker):
    """Guide to Pharmacology database linker class
    """

    def __init__(self):
        """ Initialize an object of the class GuidePharmaLinker
        """
        super(GuidePharmaLinker, self).__init__()
        self._database_name = "pharma_guide"

        self._load_linking_dictionaries()

    def convert_guide_pharmacology_to_drugbank(self, ligand_list):
        """ Convert a list of Guide to Pharmacology ligand ids to DrugBank drug ids

        Parameters
        ----------
        ligand_list : list
            a list of Guide to Pharmacology ligand ids  e.g. ['3310', '1188']

        Returns
        -------
        list
            a list of lists of DrugBank drug ids
        """
        return self._convert_ids_using_target_dictionary(ligand_list, "pharma_guide_to_drugbank")


class PubchemLinker(DatabaseLinker):
    """Pubchem database linker class
    """

    def __init__(self):
        """ Initialize an object of the class PubchemLinker
        """
        super(PubchemLinker, self).__init__()
        self._database_name = "pubchem"

        self._load_linking_dictionaries()

    def convert_substance_to_drugbank(self, substance_list):
        """ Convert a list of Pubchem subatance ids to DrugBank drug ids

        Parameters
        ----------
        substance_list : list
            a list of Pubchem subatance ids  e.g. ['46507042', '46504860']

        Returns
        -------
        list
            a list of lists of DrugBank drug ids
        """
        return self._convert_ids_using_target_dictionary(substance_list, "pubchem_substance_to_drugbank")

    def convert_substance_to_kegg_drug(self, substance_list):
        """ Convert a list of Pubchem subatance ids to KEGG drug ids

        Parameters
        ----------
        substance_list : list
            a list of Pubchem subatance ids  e.g. ['46507042', '46504860']

        Returns
        -------
        list
            a list of lists of KEGG drug ids
        """
        return self._convert_ids_using_target_dictionary(substance_list, "pubchem_substance_to_kegg_drug")

    def convert_compound_to_drugbank(self, compound_list):
        """ Convert a list of Pubchem compound ids to DrugBank drug ids

        Parameters
        ----------
        compound_list : list
            a list of Pubchem compound ids  e.g. ['45267103', '16134395']

        Returns
        -------
        list
            a list of lists of DrugBank drug ids
        """
        return self._convert_ids_using_target_dictionary(compound_list, "pubchem_compound_to_drugbank")

    def convert_substance_to_sider(self, substance_list):
        """ Convert a list of Pubchem subatance ids to SIDER drug ids

        Parameters
        ----------
        substance_list : list
            a list of Pubchem subatance ids  e.g. ['46507042', '46504860']

        Returns
        -------
        list
            a list of lists of SIDER drug ids
        """
        return self._convert_ids_using_target_dictionary(substance_list, "pubchem_substance_to_sider")

    def convert_compound_to_sider(self, compound_list):
        """ Convert a list of Pubchem compound ids to SIDER drug ids

        Parameters
        ----------
        compound_list : list
            a list of Pubchem compound ids  e.g. ['45267103', '16134395']

        Returns
        -------
        list
            a list of lists of SIDER drug ids
        """
        return self._convert_ids_using_target_dictionary(compound_list, "pubchem_compound_to_sider")


class ChebiLinker(DatabaseLinker):
    """ Chebi database linker class
    """

    def __init__(self):
        """ Initialize an object of the class ChebiLinker
        """
        super(ChebiLinker, self).__init__()
        self._database_name = "chebi"

        self._load_linking_dictionaries()


    def convert_chemical_to_drugbank_drugs(self, chemical_list):
        """ Convert a list of Chebi chemical ids to drugbank drugs

        Parameters
        ----------
        chemical_list : list
            a list of chebi chemical ids e.g. ['CHEBI:16347', 'CHEBI:11060']

        Returns
        -------
        list
            a list of lists of drugbank drugs
        """
        return self._convert_ids_using_target_dictionary(chemical_list, "chebi_to_drugbank")

    def convert_chemical_to_kegg_drugs(self, chemical_list):
        """ Convert a list of Chebi chemical ids to KEGG drugs

        Parameters
        ----------
        chemical_list : list
            a list of chebi chemical ids e.g. ['CHEBI:16347', 'CHEBI:11060']

        Returns
        -------
        list
            a list of lists of KEGG drugs
        """
        return self._convert_ids_using_target_dictionary(chemical_list, "chebi_to_kegg_drug")

    def convert_chemical_to_sider(self, chemical_list):
        """ Convert a list of Chebi chemical ids to SIDER drugs

        Parameters
        ----------
        chemical_list : list
            a list of chebi chemical ids e.g. ['CHEBI:16347', 'CHEBI:11060']

        Returns
        -------
        list
            a list of lists of SIDER compounds
        """
        return self._convert_ids_using_target_dictionary(chemical_list, "chebi_to_sider")


class DPDLinker(DatabaseLinker):
    """ Drug Product database linker class
    """

    def __init__(self):
        """ Initialize an object of the class DPDLinker
        """
        super(DPDLinker, self).__init__()
        self._database_name = "dpd"

        self._load_linking_dictionaries()


    def convert_drug_to_drugbank_drugs(self, drug_list):
        """ Convert a list of DPD drug ids to drugbank drugs

        Parameters
        ----------
        drug_list : list
            a list of DPD drug ids e.g. ['11916', '13175']

        Returns
        -------
        list
            a list of lists of drugbank drugs
        """
        return self._convert_ids_using_target_dictionary(drug_list, "dpd_to_drugbank")


class BindingDBLinker(DatabaseLinker):
    """ BindingDB database linker class
    """

    def __init__(self):
        """ Initialize an object of the class BindingDBLinker
        """
        super(BindingDBLinker, self).__init__()
        self._database_name = "bindingdb"

        self._load_linking_dictionaries()


    def convert_compound_to_drugbank_drugs(self, compound_list):
        """ Convert a list of BindingDB compound ids to drugbank drugs

        Parameters
        ----------
        compound_list : list
            a list of BindingDB compound ids e.g. ['50369395', '50022815']

        Returns
        -------
        list
            a list of lists of drugbank drugs
        """
        return self._convert_ids_using_target_dictionary(compound_list, "bindingdb_to_drugbank")

    def convert_compound_to_sider(self, compound_list):
        """ Convert a list of BindingDB compound ids to sider drugs

        Parameters
        ----------
        compound_list : list
            a list of BindingDB compound ids e.g. ['50369395', '50022815']

        Returns
        -------
        list
            a list of lists of sider drugs
        """
        return self._convert_ids_using_target_dictionary(compound_list, "bindingdb_to_sider")


class HSDBLinker(DatabaseLinker):
    """ HSDB database linker class
    """

    def __init__(self):
        """ Initialize an object of the class HSDBLinker
        """
        super(HSDBLinker, self).__init__()
        self._database_name = "hsdb"

        self._load_linking_dictionaries()


    def convert_substance_to_kegg_drugs(self, substance_list):
        """ Convert a list of HSDB substance ids to KEGG drugs

        Parameters
        ----------
        substance_list : list
            a list of HSDB substance ids e.g. ['Donepezil', 'Efalizumab']

        Returns
        -------
        list
            a list of lists of KEGG drugs
        """
        return self._convert_ids_using_target_dictionary(substance_list, "hsdb_to_kegg_drug")


class IupharLinker(DatabaseLinker):
    """ Iuphar database linker class
    """

    def __init__(self):
        """ Initialize an object of the class IupharLinker
        """
        super(IupharLinker, self).__init__()
        self._database_name = "iuphar"

        self._load_linking_dictionaries()


    def convert_ligand_to_drugbank_drugs(self, ligand_list):
        """ Convert a list of IUPHAR ligand ids to drugbank drugs

        Parameters
        ----------
        ligand_list : list
            a list of IUPHAR ligand ids e.g. ['2174', '1188']

        Returns
        -------
        list
            a list of lists of drugbank drugs
        """
        return self._convert_ids_using_target_dictionary(ligand_list, "iuphar_to_drugbank")


class KnapsackLinker(DatabaseLinker):
    """ Knapsack database linker class
    """

    def __init__(self):
        """ Initialize an object of the class KnapsackLinker
        """
        super(KnapsackLinker, self).__init__()
        self._database_name = "knapsack"

        self._load_linking_dictionaries()


    def convert_metabolite_to_kegg_drugs(self, metabolite_list):
        """ Convert a list of knapsack metabolite ids to KEGG drugs

        Parameters
        ----------
        metabolite_list : list
            a list of knapsack metabolite ids e.g. ['C00001433', 'C00007279']

        Returns
        -------
        list
            a list of lists of KEGG drugs
        """
        return self._convert_ids_using_target_dictionary(metabolite_list, "knapsack_to_kegg_drug")


class LingandBoxLinker(DatabaseLinker):
    """ LigandBox database linker class
    """

    def __init__(self):
        """ Initialize an object of the class LingandBoxLinker
        """
        super(LingandBoxLinker, self).__init__()
        self._database_name = "ligandbox"

        self._load_linking_dictionaries()


    def convert_ligand_to_kegg_drugs(self, ligand_list):
        """ Convert a list of LigandBox ligand ids to KEGG drugs

        Parameters
        ----------
        ligand_list : list
            a list of LigandBox ligand ids e.g. ['D01109', 'D06340']

        Returns
        -------
        list
            a list of lists of KEGG drugs
        """
        return self._convert_ids_using_target_dictionary(ligand_list, "ligandbox_to_kegg_drug")


class MassBankLinker(DatabaseLinker):
    """ MassBank database linker class
    """

    def __init__(self):
        """ Initialize an object of the class MassBankLinker
        """
        super(MassBankLinker, self).__init__()
        self._database_name = "massbank"

        self._load_linking_dictionaries()


    def convert_accession_to_kegg_drugs(self, accession_list):
        """ Convert a list of MassBank accession ids to KEGG drugs

        Parameters
        ----------
        accession_list : list
            a list of MassBank accession ids e.g. ['JP005783', 'WA002017']

        Returns
        -------
        list
            a list of lists of KEGG drugs
        """
        return self._convert_ids_using_target_dictionary(accession_list, "massbank_to_kegg_drug")


class NikkajiLinker(DatabaseLinker):
    """ Nikkaji database linker class
    """

    def __init__(self):
        """ Initialize an object of the class NikkajiLinker
        """
        super(NikkajiLinker, self).__init__()
        self._database_name = "nikkaji"

        self._load_linking_dictionaries()


    def convert_chemical_to_kegg_drugs(self, chemical_list):
        """ Convert a list of Nikkaji chemical ids to KEGG drugs

        Parameters
        ----------
        chemical_list : list
            a list of Nikkaji chemical ids e.g. ['J10.483C', 'J300.856H']

        Returns
        -------
        list
            a list of lists of KEGG drugs
        """
        return self._convert_ids_using_target_dictionary(chemical_list, "nikkaji_to_kegg_drug")


class PDBLinker(DatabaseLinker):
    """ PDB database linker class
    """

    def __init__(self):
        """ Initialize an object of the class PDBLinker
        """
        super(PDBLinker, self).__init__()
        self._database_name = "pdb"

        self._load_linking_dictionaries()

    def convert_chemical_to_drugbank_drugs(self, chemical_list):
        """ Convert a list of PDB chemical ids to Drugbank drugs

        Parameters
        ----------
        chemical_list : list
            a list of PDB chemical ids e.g. ['PLP', 'HIS']

        Returns
        -------
        list
            a list of lists of Drugbank drugs
        """
        return self._convert_ids_using_target_dictionary(chemical_list, "pdb_to_drugbank")

    def convert_chemical_to_kegg_drugs(self, chemical_list):
        """ Convert a list of PDB chemical ids to KEGG drugs

        Parameters
        ----------
        chemical_list : list
            a list of PDB chemical ids e.g. ['PLP', 'HIS']

        Returns
        -------
        list
            a list of lists of KEGG drugs
        """
        return self._convert_ids_using_target_dictionary(chemical_list, "pdb_to_kegg_drug")


class TTDLinker(DatabaseLinker):
    """ TTD database linker class
    """

    def __init__(self):
        """ Initialize an object of the class TTDLinker
        """
        super(TTDLinker, self).__init__()
        self._database_name = "ttd"

        self._load_linking_dictionaries()

    def convert_drud_to_drugbank_drugs(self, drug_list):
        """ Convert a list of TTD drug ids to Drugbank drugs

        Parameters
        ----------
        drug_list : list
            a list of TTD drug ids e.g. ['DNC000788', 'DAP000020']

        Returns
        -------
        list
            a list of lists of Drugbank drugs
        """
        return self._convert_ids_using_target_dictionary(drug_list, "ttd_to_drugbank")


class WikipediaLinker(DatabaseLinker):
    """ Wikipedia database linker class
    """

    def __init__(self):
        """ Initialize an object of the class WikipediaLinker
        """
        super(WikipediaLinker, self).__init__()
        self._database_name = "wikipedia"

        self._load_linking_dictionaries()

    def convert_entry_to_drugbank_drugs(self, entry_list):
        """ Convert a list of Wikipedia entry ids to Drugbank drugs

        Parameters
        ----------
        entry_list : list
            a list of Wikipedia entry ids e.g. ['Dornase_alfa', 'Leuprolide']

        Returns
        -------
        list
            a list of lists of Drugbank drugs
        """
        return self._convert_ids_using_target_dictionary(entry_list, "wikipedia_to_drugbank")


class CellosaurusLinker(DatabaseLinker):
    """ Cellosaurus database linker class
    """

    def __init__(self):
        """ Initialize an object of the class CellosaurusLinker
        """
        super(CellosaurusLinker, self).__init__()
        self._database_name = "cellosaurus"

        self._load_linking_dictionaries()

    def convert_cell_line_to_name(self, cellline_list):
        """ Convert a list of Cellosaurus cell line ids to their names

        Parameters
        ----------
        cellline_list : list
            a list of Cellosaurus cell ids e.g. ['CVCL_U602', 'CVCL_0023']

        Returns
        -------
        list
            a list of lists of names
        """
        return self._convert_ids_using_target_dictionary(cellline_list, "cellosaurus_to_names")

    @property
    def cell_line_ids(self):
        """
        All available cellline ids in the database
        """
        return self._get_keys_of_target_dictionary("cellosaurus_to_names")


class MESHLinker(DatabaseLinker):
    """ MESH database linker class
    """

    def __init__(self):
        """ Initialize an object of the class MESHLinker
        """
        super(MESHLinker, self).__init__()
        self._database_name = "mesh"

        self._load_linking_dictionaries()

    def convert_disease_to_kegg_disease(self, disease_list):
        """ Convert a list of MESH disease ids to KEGG disease ids

        Parameters
        ----------
        disease_list : list
            a list of MESH disease ids e.g. ['D054198', 'D054218']

        Returns
        -------
        list
            a list of lists of KEGG disease ids
        """
        return self._convert_ids_using_target_dictionary(disease_list, "mesh_to_kegg_disease")

    def convert_disease_to_omim(self, disease_list):
        """ Convert a list of MESH disease ids to omim ids

        Parameters
        ----------
        disease_list : list
            a list of MESH disease ids e.g. ['D054198', 'D054218']

        Returns
        -------
        list
            a list of lists of omim ids
        """
        return self._convert_ids_using_target_dictionary(disease_list, "mesh_to_omim")

class ICD10Linker(DatabaseLinker):
    """ ICD10 database linker class
    """

    def __init__(self):
        """ Initialize an object of the class ICD10Linker
        """
        super(ICD10Linker, self).__init__()
        self._database_name = "icd10"

        self._load_linking_dictionaries()

    def convert_disease_to_kegg_disease(self, disease_list):
        """ Convert a list of MESH disease ids to KEGG disease ids

        Parameters
        ----------
        disease_list : list
            a list of MESH disease ids e.g. ['C91.0', 'C83.5']

        Returns
        -------
        list
            a list of lists of KEGG disease ids
        """
        return self._convert_ids_using_target_dictionary(disease_list, "icd_10_to_kegg_disease")


class ICD11Linker(DatabaseLinker):
    """ ICD11 database linker class
    """

    def __init__(self):
        """ Initialize an object of the class ICD11Linker
        """
        super(ICD11Linker, self).__init__()
        self._database_name = "icd11"

        self._load_linking_dictionaries()

    def convert_disease_to_kegg_disease(self, disease_list):
        """ Convert a list of MESH disease ids to KEGG disease ids

        Parameters
        ----------
        disease_list : list
            a list of MESH disease ids e.g. ['2A70', '2A71']

        Returns
        -------
        list
            a list of lists of KEGG disease ids
        """
        return self._convert_ids_using_target_dictionary(disease_list, "icd_11_to_kegg_disease")


