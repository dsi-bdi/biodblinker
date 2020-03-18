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

    @abstractmethod
    def _load_linking_dictionaries(self):
        """ Load database linking dictionaries
        """
        raise NotImplementedError("Method is not implemented")

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
                outcome_ids.append("")
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

    def _load_linking_dictionaries(self):
        """ Load database linking dictionaries
        """
        self._database_linking_dictionaries_paths = get_database_mapping_sources(self._database_name)

        for mapping_name, mapping_filepath in self._database_linking_dictionaries_paths.items():
            self._database_linking_dictionaries[mapping_name] = load_data_map(mapping_filepath)

    def convert_geneid_to_uniprot(self, geneid_list):
        """ Convert a list of KEGG gene ids to uniprot accessions

        Parameters
        ----------
        geneid_list : list
            a list of KEGG gene ids e.g. ['hsa:8209', 'hsa:1564']

        Returns
        -------
        list
            a list of uniprot accession codes
        """
        return self._convert_ids_using_target_dictionary(geneid_list, "geneid_to_uniprot")

    def convert_geneid_to_names(self, geneid_list):
        """ Convert a list of KEGG gene ids to gene names

        Parameters
        ----------
        geneid_list : list
            a list of KEGG gene ids e.g. ['hsa:8209', 'hsa:1564']

        Returns
        -------
        list
            a list of gene names
        """
        return self._convert_ids_using_target_dictionary(geneid_list, "geneid_to_name")

    @property
    def gene_ids(self):
        """
        All available gene ids in the database
        """
        return self._get_keys_of_target_dictionary("geneid_to_name")

    @property
    def drug_ids(self):
        """
        All available drug ids in the database
        """
        return self._get_keys_of_target_dictionary("geneid_to_name")
