from os import makedirs, remove, walk
from os.path import join, isdir, exists, isfile
from biodblinker.fileio import *
from biodblinker.config import file_open, get_data_directory, get_all_mappings_sources
from tqdm import tqdm
import requests
from collections import defaultdict
from zipfile import ZipFile
import xml.etree.ElementTree as ET
import gzip
import shutil
import time



class Species:

    def __init__(self, code, kegg_organism, node, scientific_name):
        self.code = code
        self.kegg_organism = kegg_organism
        self.node = node
        self.scientific_name = scientific_name


VALID_SPECIES = [
    Species('ARATH', 'ath', '3702', 'arabidopsis thaliana'),
    Species('BACSU', 'bsu', '224308', 'bacillus subtilis (strain 168)'),
    Species('BOVIN', 'bta', '9913', 'bos taurus'),
    Species('CAEEL', 'cel', '6239', 'caenorhabditis elegans'),
    Species('CHICK', 'gga', '9031', 'gallus gallus'),
    Species('DANRE', 'dre', '7955', 'danio rerio'),
    Species('DICDI', 'ddi', '44689', 'dictyostelium discoideum'),
    Species('DROME', 'dme', '7227', 'drosophila melanogaster'),
    Species('ECO57', 'ece', '83334', 'escherichia coli o157:h7'),
    Species('ECOLI', 'eco', '83333', 'escherichia coli (strain k12)'),
    Species('HUMAN', 'hsa', '9606', 'homo sapiens'),
    Species('MOUSE', 'mmu', '10090', 'mus musculus'),
    Species('MYCTO', 'mtc', '83331', 'mycobacterium tuberculosis (strain cdc 1551 / oshkosh)'),
    Species('MYCTU', 'mtu', '83332', 'mycobacterium tuberculosis (strain atcc 25618 / h37rv)'),
    Species('ORYSJ', 'osa', '39947', 'oryza sativa subsp. japonica'),
    Species('PONAB', 'pon', '9601', 'pongo abelii'),
    Species('RAT', 'rno', '10116', 'rattus norvegicus'),
    Species('SCHPO', 'spo', '284812', 'schizosaccharomyces pombe (strain 972 / atcc 24843)'),
    Species('XENLA', 'xla', '8355', 'xenopus laevis'),
    Species('YEAST', 'sce', '559292', 'saccharomyces cerevisiae (strain atcc 204508 / s288c)'),
    Species('PIG', 'ssc', '9823', 'Sus scrofa')
]

class MappingGenerator():

    def __init__(self):
        """
        """
        self._data_dir = join(get_data_directory(), 'data')
        self._source_dir = join(get_data_directory(), 'sources')

    def _download_uniprot_sources(self, sources_dir):
        """ Download uniprot mappings file

        Parameters
        ----------
        sources_dir : str
            the path to output the mappings

        Returns
        --------
        str
            the path to the mappings file
        """
        swissprot_file_url = 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz'
        uniprot_mapping_file_url = "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz"

        uniprot_mapping_filepath = join(sources_dir, "./uniprot_mappings.dat.gz")
        swissprot_mapping_filepath = join(sources_dir, "./swissprot.xml.gz")

        download_file_md5_check(uniprot_mapping_file_url, uniprot_mapping_filepath)
        download_file_md5_check(swissprot_file_url, swissprot_mapping_filepath)
        return uniprot_mapping_filepath, swissprot_mapping_filepath

    def _get_included_accs(self, swissprot_file):
        """ Get the set of swissprot accessions to include

        Parameters
        ----------
        swissprot_file : str
            the path to the swissprot xml file

        Returns
        --------
        set
            The set of swissprot accessions to include
        """

        ns = {'up': 'http://uniprot.org/uniprot'}
        valid_species = list(map(lambda x: x.code, VALID_SPECIES))
        with file_open(swissprot_file) as swissprot_fd:
            valid_accs = set()
            for _, elem in tqdm(ET.iterparse(swissprot_fd), 'Processing swissprot mapping'):
                if elem.tag == '{http://uniprot.org/uniprot}entry':
                    uniprot_acc = elem.find('./up:accession', ns).text
                    name = elem.find('./up:name', ns).text

                    species = name.split('_')[-1]
                    valid_accs.add(uniprot_acc)
                    # if species in valid_species:
                    #     valid_accs.add(uniprot_acc)

                    elem.clear()
        return valid_accs

    def _map_uniprot(self, uniprot_mapping_file, valid_accs):
        """ Map uniprot proteins to several databases

        Parameters
        ----------
        uniprot_mapping_file : str
            the path to the uniprot mappings file

        Returns
        --------
        dict
            dictionary of dictionaries mapping uniprot proteins to target databases 
        dict
            dictionary of dictionaries mapping target databases to uniprot proteins
        """
        target_databases = {"Gene_Name", "PharmGKB", "STRING",
                            "HGNC", "BioGrid", "RefSeq", "UniParc", "KEGG",
                            "Ensembl", "MIM", "MINT", "UniProtKB-ID", "EMBL-CDS",
                            "DIP", "ProteomicsDB", "UniRef100", "GI", "GeneDB",
                            "EMBL", "ChEMBL", "KO", "MIM", "Gene_Synonym"}

        db_mapping_dictionaries_dict = {dbname: dict() for dbname in target_databases}
        db_mapping_dictionaries_reverse_dict = {dbname: dict() for dbname in target_databases}
        with file_open(uniprot_mapping_file) as uniprot_map_fd:
            for line in tqdm(uniprot_map_fd, desc="Processing UNIPROT mapping file."):
                uniprot_acc, db_tag, db_id = line.strip().split("\t")
                if uniprot_acc not in valid_accs:
                    continue
                if db_tag in target_databases:
                    if uniprot_acc not in db_mapping_dictionaries_dict[db_tag]:
                        db_mapping_dictionaries_dict[db_tag][uniprot_acc] = [db_id]
                    else:
                        db_mapping_dictionaries_dict[db_tag][uniprot_acc].append(db_id)

                    if db_id not in db_mapping_dictionaries_reverse_dict[db_tag]:
                        db_mapping_dictionaries_reverse_dict[db_tag][db_id] = [uniprot_acc]
                    else:
                        db_mapping_dictionaries_reverse_dict[db_tag][db_id].append(uniprot_acc)

        return db_mapping_dictionaries_dict, db_mapping_dictionaries_reverse_dict

    def _export_uniprot(self, db_mapping_dictionaries_dict, uniprot_mappings_dp):
        """ Exprot the uniprot mappings

        Parameters
        ----------
        db_mapping_dictionaries_dict : dict
            dictionary of dictionaries mapping uniprot proteins to target databases 
        uniprot_mappings_dp : str
            the directory to output the mappings
        """
        for db_name in tqdm(db_mapping_dictionaries_dict, desc="Exporting uniprot target databases mappings"):
            db_map_fp = join(uniprot_mappings_dp, f"acc_to_{db_name.lower()}.txt")
            db_map_fd = open(db_map_fp, "w", encoding='utf-8')
            for key, val in db_mapping_dictionaries_dict[db_name].items():
                val_txt = ";".join(val)
                db_map_fd.write(f"{key}\t{val_txt}\n")
            db_map_fd.close()

    def _export_uniprot_reverse(self, db_mapping_dictionaries_reverse_dict, data_dir):
        """ Exprot the reversed uniprot mappings

        Parameters
        ----------
        db_mapping_dictionaries_reverse_dict : dict
            dictionary of dictionaries mapping target databases to uniprot proteins
        data_dir : str
            the directory to output the mappings
        """
        # change "Gene_Synonym" to "Gene_Name" to treat both equally
        # add gene synonyms to gene name database entry in the reverse maps dictionary and remove the synonyms entry
        db_mapping_dictionaries_reverse_dict["Gene_Name"].update(db_mapping_dictionaries_reverse_dict["Gene_Synonym"].copy())
        del db_mapping_dictionaries_reverse_dict["Gene_Synonym"]

        # export to files in respective dictionaries
        for db_name in tqdm(db_mapping_dictionaries_reverse_dict, desc="Exporting uniprot target databases reverse mappings"):
            db_mappings_dp = join(data_dir, db_name.lower())
            makedirs(db_mappings_dp) if not isdir(db_mappings_dp) else None
            db_rev_map_fp = join(db_mappings_dp, f"{db_name.lower()}_to_uniprot.txt")
            db_rev_map_fd = open(db_rev_map_fp, "w", encoding='utf-8')
            for key, val in db_mapping_dictionaries_reverse_dict[db_name].items():
                val_txt = ";".join(val)
                db_rev_map_fd.write(f"{key}\t{val_txt}\n")
            db_rev_map_fd.close()

    def generate_uniprot_mappings(self):
        """ Generate mappings from uniprot proteins """

        sources_dp = self._source_dir
        mappings_dp = join(self._data_dir, 'uniprot')

        makedirs(sources_dp) if not isdir(sources_dp) else None
        makedirs(mappings_dp) if not isdir(mappings_dp) else None

        uniprot_mapping_file, swissprot_file = self._download_uniprot_sources(sources_dp)
        uniprot_accs = self._get_included_accs(swissprot_file)

        uniprot_mappings, uniprot_mappings_rev = self._map_uniprot(uniprot_mapping_file, uniprot_accs)
        self._export_uniprot(uniprot_mappings, mappings_dp)
        self._export_uniprot_reverse(uniprot_mappings_rev, self._data_dir)

        with open(join(mappings_dp, 'uniprot_names.txt'), 'w', encoding='utf-8') as writer:
            for acc, ids in uniprot_mappings['UniProtKB-ID'].items():
                names = set()
                for pid in ids:
                    name = "_".join(pid.split('_')[:-1])
                    names.add(name)
                writer.write(f'{acc}\t{";".join(names)}\n')

    def _map_kegg_names(self, database, mappings_dp, file_mode='w', is_gene_db=False, database_name=None):
        """ map kegg entities to their names

        Parameters
        ----------
        database : str
            kegg entity database
        mappings_dp : str
            the directory to output the name mappings

        Returns
        -------
        set
            the set of entity ids in the specified database
        """
        id_set = set()

        if database_name is None:
            database_name = database

        resp = requests.get(f'http://rest.kegg.jp/list/{database}')

        if resp.ok:
            with open(join(mappings_dp, f'{database_name}_names.txt'), file_mode, encoding='utf-8') as writer:
                for line in resp.iter_lines(decode_unicode=True):
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        (source_str, name_str) = parts[0:2]
                        _, did = source_str.split(':')
                        if did[0].isdigit():
                            did = source_str
                        elif is_gene_db:
                            did = source_str
                        id_set.add(did)

                        names = set()
                        for name in name_str.split(';')[0].split(','):
                            name = name.strip()
                            if database == 'drug' or database == 'environw':
                                index = name.rfind('(')
                                if index > -1:
                                    name = name[:index].strip()

                            names.add(name)
                        writer.write(f'{did}\t{";".join(names)}\n')
                    else:
                        source_str = parts[0]
                        _, did = source_str.split(':')
                        if did[0].isdigit():
                            did = source_str
                        id_set.add(did)
                        writer.write(f'{did}\t-\n')

        return id_set

    def _link_kegg_db(self, database, target):
        """ link kegg entities to specified target

        Parameters
        ----------
        database : str
            kegg entity database
        target : str
            database to link to

        Returns
        -------
        dict
            A dictionary mapping entities in the database to the target database
        """
        mapping_dict = defaultdict(set)
        resp = requests.get(f'http://rest.genome.jp/link/{target}/{database}')
        if resp.ok:
            for line in resp.iter_lines(decode_unicode=True):
                source_str, target_str, _ = line.strip().split('\t')
                _, did = source_str.split(':')
                if did[0].isdigit():
                    did = source_str

                _, tid = target_str.split(':')
                if target == 'sider':
                    sider_base = 'CID100000000'
                    tid = sider_base[:len(sider_base)-len(tid)]+tid
                mapping_dict[did].add(tid)

        return mapping_dict

    def _download_kegg_diseases(self, sources_dp):
        """ Download kegg medicus diseases

        Parameters
        ----------
        sources_dp : str
            the path to output the diseases file

        Returns
        --------
        str
            the path to the diseases file
        """
        diseases_file_url = "ftp://ftp.genome.jp/pub/kegg/medicus/disease/disease"
        diseases_filepath = join(sources_dp, "./kegg_diseases.txt")
        download_file_md5_check(diseases_file_url, diseases_filepath)

        return diseases_filepath

    def parse_disease_mappings(self, diseases_fp, mappings_dp):
        disease_mapping_dict = {'mesh': {}, 'icd_10': {}, 'icd_11': {}}
        rev_disease_mapping_dict = {'mesh': {}, 'icd_10': {}, 'icd_11': {}}
        disease_id = ''
        current_section = ''
        with open(diseases_fp, 'r') as diseases_fd:
            for line in diseases_fd:
                section = line.split(' ')[0]
                # New Section
                if len(section.strip()) > 0:
                    current_section = section.strip()

                if line.startswith('///'):
                    disease_id = ''
                elif current_section == 'ENTRY':
                    disease_id = line.split()[1]
                elif current_section == 'DBLINKS':
                    if len(disease_id) == 0:
                        raise Exception('DiseseID not set for link')

                    line = line.replace(section, '')
                    parts = line.split()
                    db = parts[0]
                    target_ids = parts[1:]

                    if db == 'ICD-11:':
                        if disease_id not in disease_mapping_dict['icd_11']:
                            disease_mapping_dict['icd_11'][disease_id] = []
                        disease_mapping_dict['icd_11'][disease_id].extend(target_ids)
                        for target in target_ids:
                            if target not in rev_disease_mapping_dict['icd_11']:
                                rev_disease_mapping_dict['icd_11'][target] = []
                            rev_disease_mapping_dict['icd_11'][target].append(disease_id)
                    elif db == 'ICD-10:':
                        if disease_id not in disease_mapping_dict['icd_10']:
                            disease_mapping_dict['icd_10'][disease_id] = []
                        disease_mapping_dict['icd_10'][disease_id].extend(target_ids)
                        for target in target_ids:
                            if target not in rev_disease_mapping_dict['icd_10']:
                                rev_disease_mapping_dict['icd_10'][target] = []
                            rev_disease_mapping_dict['icd_10'][target].append(disease_id)
                    elif db == 'MeSH:':
                        if disease_id not in disease_mapping_dict['mesh']:
                            disease_mapping_dict['mesh'][disease_id] = []

                        for target in target_ids:
                            disease_mapping_dict['mesh'][disease_id].append(target)
                            if target not in rev_disease_mapping_dict['mesh']:
                                rev_disease_mapping_dict['mesh'][target] = []
                            rev_disease_mapping_dict['mesh'][target].append(disease_id)

            for db, mapping in disease_mapping_dict.items():
                with open(join(mappings_dp, f'disease_to_{db}.txt'), 'w', encoding='utf-8') as writer:
                    for disease, map_ids in mapping.items():
                        map_ids = set(map_ids)
                        if len(map_ids) == 0:
                            map_ids.add('-')
                        else:
                            writer.write(f'{disease}\t{";".join(map_ids)}\n')

            for db, mapping in rev_disease_mapping_dict.items():
                target_dp = join(self._data_dir, db)
                makedirs(target_dp) if not isdir(target_dp) else None
                with open(join(target_dp, f'{db}_to_kegg_disease.txt'), 'w', encoding='utf-8') as writer:
                    for target_id, map_ids in mapping.items():
                        map_ids = set(map_ids)
                        if len(map_ids) == 0:
                            map_ids.add('-')
                        else:
                            writer.write(f'{target_id}\t{";".join(map_ids)}\n')

    def generate_kegg_mappings(self):
        """ Generate mappings for kegg entities """
        gene_xref = ['uniprot', 'ensembl']

        species_list = map(lambda x: x.kegg_organism, VALID_SPECIES)
        kegg_source_targets = {
            'drug':
            [
                'hsdb', 'hmdb', 'nikkaji', 'chembl',
                'knapsack', 'pubchem', 'chebi', 'pdb-ccd',
                'lipidbank', 'lipidmaps', 'ligandbox', 'massbank'
            ],
            'disease':
            [
                'omim', 'pubmed'
            ],
            'pathway': [],
            'glycan': [],
            'network': []
        }

        sources_dp = self._source_dir
        mappings_dp = join(self._data_dir, 'kegg')
        makedirs(sources_dp) if not isdir(sources_dp) else None
        makedirs(mappings_dp) if not isdir(mappings_dp) else None

        diseases_fp = self._download_kegg_diseases(sources_dp)
        self.parse_disease_mappings(diseases_fp, mappings_dp)

        for source, targets in tqdm(kegg_source_targets.items(), 'Processing KEGG mappings'):
            source_ids = self._map_kegg_names(source, mappings_dp)
            # Sleep between requests
            time.sleep(0.2)
            for target in targets:
                mapping_dict = self._link_kegg_db(source, target)
                # Sleep between requests
                time.sleep(0.2)
                if len(mapping_dict) > 0:
                    source_name = source

                    reverse_map = defaultdict(set)
                    with open(join(mappings_dp, f'{source_name}_to_{target}.txt'), 'w', encoding='utf-8') as writer:
                        for sid in source_ids:
                            map_ids = mapping_dict[sid]
                            if len(map_ids) == 0:
                                map_ids.add('-')
                            else:
                                for mapped_id in map_ids:
                                    reverse_map[mapped_id].add(sid)
                            writer.write(f'{sid}\t{";".join(map_ids)}\n')

                    target_dp = join(self._data_dir, target)
                    makedirs(target_dp) if not isdir(target_dp) else None
                    with open(join(target_dp, f'{target}_to_kegg_{source_name}.txt'), 'w', encoding='utf-8') as writer:
                        for tid, kids in reverse_map.items():
                            writer.write(f'{tid}\t{";".join(kids)}\n')

        file_mode = 'w'
        for species in tqdm(species_list, 'Processing KEGG gene mappings'):
            source_ids = self._map_kegg_names(species, mappings_dp, file_mode, True, database_name='gene')
            # Sleep between requests
            time.sleep(0.2)
            for target in gene_xref:
                mapping_dict = self._link_kegg_db(species, target)
                # Sleep between requests
                time.sleep(0.2)
                if len(mapping_dict) > 0:
                    source_name = 'gene'

                    reverse_map = defaultdict(set)
                    with open(join(mappings_dp, f'{source_name}_to_{target}.txt'), file_mode, encoding='utf-8') as writer:
                        for sid in source_ids:
                            map_ids = mapping_dict[sid]
                            if len(map_ids) == 0:
                                map_ids.add('-')
                            else:
                                for mapped_id in map_ids:
                                    reverse_map[mapped_id].add(sid)
                            writer.write(f'{sid}\t{";".join(map_ids)}\n')

                    target_dp = join(self._data_dir, target)
                    makedirs(target_dp) if not isdir(target_dp) else None
                    with open(join(target_dp, f'{target}_to_kegg_{source_name}.txt'), file_mode, encoding='utf-8') as writer:
                        for tid, kids in reverse_map.items():
                            writer.write(f'{tid}\t{";".join(kids)}\n')
                # After the first iteration append to the files instead of overwriting
                file_mode = 'a'

        kegg_to_mesh = {}
        kegg_to_omim = {}
        with open(join(mappings_dp, 'disease_to_mesh.txt'), 'r') as fd:
            for line in fd:
                if line.strip().split('\t') == '-':
                    continue
                kegg_id, mesh_ids = line.strip().split('\t')
                mesh_ids = mesh_ids.split(';')
                kegg_to_mesh[kegg_id] = mesh_ids
        with open(join(mappings_dp, 'disease_to_omim.txt'), 'r') as fd:
            for line in fd:
                if line.strip().split('\t') == '-':
                    continue
                kegg_id, omim_ids = line.strip().split('\t')
                omim_ids = omim_ids.split(';')
                kegg_to_omim[kegg_id] = omim_ids

        with open(join(self._data_dir, 'mesh', 'mesh_to_omim.txt'), 'w') as fd:
            for kegg_id, mesh_mappings in kegg_to_mesh.items():
                if kegg_id in kegg_to_omim:
                    for mesh_id in mesh_mappings:
                        fd.write(f'{mesh_id}\t{";".join(kegg_to_omim[kegg_id])}\n')

        with open(join(self._data_dir, 'omim', 'omim_to_mesh.txt'), 'w') as fd:
            for kegg_id, omim_mappings in kegg_to_omim.items():
                if kegg_id in kegg_to_mesh:
                    for omim_id in omim_mappings:
                        fd.write(f'{omim_id}\t{";".join(kegg_to_mesh[kegg_id])}\n')

    def _download_drugbank(self, sources_dp, username, password):
        """ Download drugbank

        Parameters
        ----------
        sources_dp : str
            the path to output the drugbank file
        username : str
            drugbank username
        password : str
            drugbank password

        Returns
        --------
        str
            the path to the drugbank file
        """
        drugbank_file_url = "https://www.drugbank.ca/releases/latest/downloads/all-full-database"
        drugbank_filepath = join(sources_dp, "./drugbank_all_full_database.xml.zip")
        download_file_md5_check(drugbank_file_url, drugbank_filepath, username=username, password=password)

        return drugbank_filepath

    def _map_drugbank(self, drugbank_filepath, drugbank_targets):
        """ map drugbank drugs to the target databases

        Parameters
        ----------
        drugbank_filepath : str
            the path to the drugbank file
        drugbank_targets : dict
            dictionary of target databases mapped to normalized names

        Returns
        -------
        set
            the set of drugbank drug ids
        dict
            dictionary of dictionaries mapping drugbank drugs to the target databases
        """
        id_set = set()
        drugbank_targets_map = {}
        drugbank_targets_map['names'] = defaultdict(set)
        drugbank_targets_map['ids'] = defaultdict(set)
        for target in drugbank_targets.values():
            drugbank_targets_map[target] = defaultdict(set)

        ns = {'db': 'http://www.drugbank.ca'}
        with ZipFile(drugbank_filepath, 'r') as dbzip:
            with dbzip.open('full database.xml', force_zip64=True) as xmlfile:
                for _, elem in tqdm(ET.iterparse(xmlfile), 'Processing Drugbank mapping'):
                    # Check the length of the drug element as pathways also contain drug elements
                    if elem.tag == '{http://www.drugbank.ca}drug' and len(elem) > 2:
                        drug_id_elem = elem.find('./db:drugbank-id[@primary="true"]', ns)
                        did = drug_id_elem.text

                        id_set.add(did)
                        name_elem = elem.find('./db:name', ns)

                        if name_elem is not None:
                            drugbank_targets_map['names'][did].add(name_elem.text.strip())

                        synonyms = elem.find('./db:synonyms', ns)
                        if synonyms is not None:
                            for synonym in synonyms:
                                drugbank_targets_map['names'][did].add(synonym.text.strip())

                        id_elems = elem.findall('./db:drugbank-id', ns)

                        if id_elems is not None:
                            for id in id_elems:
                                if id != did:
                                    drugbank_targets_map['ids'][did].add(id.text)

                        ext_ids = elem.findall('./db:external-identifiers/db:external-identifier', ns)
                        for ext_id in ext_ids:
                            source = ext_id.find('./db:resource', ns)
                            if source.text in drugbank_targets:
                                target = drugbank_targets[source.text]
                                target_id = ext_id.find('./db:identifier', ns).text
                                drugbank_targets_map[target][did].add(target_id)

                        elem.clear()

        return id_set, drugbank_targets_map

    def generate_drugbank_mappings(self, username, password):
        """ Generate mappings for drugbank drugs

        Parameters
        ----------
        username : str
            drugbank username
        password : str
            drugbank password
        """
        drugbank_targets = {
            'BindingDB': 'bindingdb',
            'ChEBI': 'chebi',
            'ChEMBL': 'chembl',
            'ChemSpider': 'chemspider',
            'Drugs Product Database (DPD)': 'dpd',
            'IUPHAR': 'iuphar',
            'KEGG Drug': 'kegg_drug',
            'PDB': 'pdb',
            'PDRhealth': 'pdrhealth',
            'PharmGKB': 'pharmgkb',
            'PubChem Compound': 'pubchem_compound',
            'PubChem Substance': 'pubchem_substance',
            'Therapeutic Targets Database': 'ttd',
            'Guide to Pharmacology': 'pharma_guide',
            'UniProtKB': 'uniprot',
            'Wikipedia': 'wikipedia',
            'GenBank': 'genbank'
        }

        sources_dp = self._source_dir
        mappings_dp = join(self._data_dir, 'drugbank')

        makedirs(sources_dp) if not isdir(sources_dp) else None
        makedirs(mappings_dp) if not isdir(mappings_dp) else None

        drugbank_filepath = self._download_drugbank(sources_dp, username, password)
        drugbank_ids, drugbank_targets_map = self._map_drugbank(drugbank_filepath, drugbank_targets)
        for target, target_map in tqdm(drugbank_targets_map.items(), 'Exporting Drugbank mappings'):
            if len(target_map) == 0:
                continue
            reverse_map = defaultdict(set)
            with open(join(mappings_dp, f'drugbank_to_{target}.txt'), 'w', encoding='utf-8') as writer:
                for did in sorted(drugbank_ids):
                    mapped = target_map[did]
                    if len(mapped) == 0:
                        mapped.add('-')
                    else:
                        for mapped_id in mapped:
                            reverse_map[mapped_id].add(did)
                    writer.write(f'{did}\t{";".join(mapped)}\n')

            target_dp = join(self._data_dir, target)
            makedirs(target_dp) if not isdir(target_dp) else None

            with open(join(target_dp, f'{target}_to_drugbank.txt'), 'w', encoding='utf-8') as writer:
                for tid, dids in reverse_map.items():
                    writer.write(f'{tid}\t{";".join(dids)}\n')

    def _download_sider_files(self, source_dir):
        """ Download sider names file and stitch mappings file

        Parameters
        ----------
        sources_dir : str
            the path to save the files

        Returns
        -------
        str
            the path to the sider names file
        str
            the path to the stitch mappings file
        """
        sider_names_file = 'http://sideeffects.embl.de/media/download/drug_names.tsv'
        sider_names_filepath = join(source_dir, "./sider_drug_names.tsv")
        download_file_md5_check(sider_names_file, sider_names_filepath)
        stitch_mapping_file = 'http://stitch.embl.de/download/chemical.sources.v5.0.tsv.gz'
        stitch_mapping_filepath = join(source_dir, "./stitch_mapping.tsv.gz")
        download_file_md5_check(stitch_mapping_file, stitch_mapping_filepath)
        return sider_names_filepath, stitch_mapping_filepath

    def _map_sider(self, sider_targets, names_fp, mappings_fp):
        """ Map sider drugs to their names and the specified target databases

        Parameters
        ----------
        sider_targets : dict
            dictionary mapping target databses to their normalized names
        names_fp : str
            the path to the sider names file
        mappings_fp : str
            the path to the stitch mappings file

        Returns
        -------
        set
            the set of sider drug ids
        dict
            dictionary of dictionaries mapping sider drugs to the target databases
        """
        sider_ids = set()
        sider_target_map = {}
        sider_target_map['names'] = defaultdict(set)
        for _, target in sider_targets.items():
            sider_target_map[target] = defaultdict(set)

        with open(names_fp, 'r') as names_fd:
            for line in names_fd:
                sid, name = line.strip().split('\t')
                sider_target_map['names'][sid].add(name)
                sider_ids.add(sid)

        with gzip.open(mappings_fp, 'rt') as stitch_fd:
            next(stitch_fd)
            for line in tqdm(stitch_fd, 'Processing Sider mapping'):
                if line.startswith('#'):
                    continue
                (mole, _, target, target_id) = line.strip().split('\t')
                sid = 'CID1'+mole[4:]
                if target in sider_targets and sid in sider_ids:
                    target_name = sider_targets[target]
                    sider_target_map[target_name][sid].add(target_id)

        return sider_ids, sider_target_map

    def generate_sider_mappings(self):
        """ Generate mappings for sider drugs"""
        sider_targets = {
            'BindingDB': 'bindingdb',
            'ChEBI': 'chebi',
            'ChEMBL': 'chembl',
            'DrugBank': 'drugbank',
            'KEGG': 'kegg',
            'PC': 'pubchem_compound',
            'PS': 'pubchem_substance',
            'ATC': 'atc'
        }

        sources_dp = self._source_dir
        mappings_dp = join(self._data_dir, 'sider')

        makedirs(sources_dp) if not isdir(sources_dp) else None
        makedirs(mappings_dp) if not isdir(mappings_dp) else None

        names_fp, mappings_fp = self._download_sider_files(sources_dp)
        sider_ids, sider_target_map = self._map_sider(sider_targets, names_fp, mappings_fp)

        for target, target_map in tqdm(sider_target_map.items(), 'Exporting Sider mapping'):
            if len(target_map) == 0:
                continue

            reverse_map = defaultdict(set)
            with open(join(mappings_dp, f'sider_to_{target}.txt'), 'w', encoding='utf-8') as writer:
                for sid in sorted(sider_ids):
                    mapped = target_map[sid]
                    if len(mapped) == 0:
                        mapped.add('-')
                    else:
                        for mapped_id in mapped:
                            reverse_map[mapped_id].add(sid)
                    writer.write(f'{sid}\t{";".join(mapped)}\n')

            target_dp = join(self._data_dir, target)
            makedirs(target_dp) if not isdir(target_dp) else None

            with open(join(target_dp, f'{target}_to_sider.txt'), 'w', encoding='utf-8') as writer:
                for tid, sids in reverse_map.items():
                    writer.write(f'{tid}\t{";".join(sids)}\n')

    def _download_cellosaurus_files(self, source_dir):
        """ Download cellosaurus file

        Parameters
        ----------
        sources_dir : str
            the path to save the files

        Returns
        -------
        str
            the path to the cellosaurus file
        """
        cello_file_url = 'ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml'
        cello_filepath = join(source_dir, "./cellosaurus.xml")
        download_file_md5_check(cello_file_url, cello_filepath)
        return cello_filepath

    def _map_cello_names(self, cello_fp, mapping_fp):
        """ Map cellosaurus accessions to their names

        Parameters
        ----------
        cello_fp : str
            the path to the cellosaurus file
        mapping_fp : str
            the path to the output file
        """
        with open(mapping_fp, 'w') as writer:
            with open(cello_fp) as xmlfile:
                for _, elem in tqdm(ET.iterparse(xmlfile), 'Processing Cellosaurus mapping'):
                    if elem.tag == 'cell-line' and len(elem) > 2:
                        names = elem.findall('./name-list/name')
                        namelist = []
                        if names is not None:
                            for name in names:
                                namelist.append(name.text.strip())
                        if len(names) == 0:
                            continue
                        name_str = ";".join(namelist)
                        accessions = elem.findall('./accession-list/accession')
                        if accessions is not None:
                            for accession in accessions:
                                writer.write(f'{accession.text.strip()}\t{name_str}\n')
                        elem.clear()

    def generate_cellosaurus_mappings(self):
        """ Generate mappings for cellosaurus celllines"""
        sources_dp = self._source_dir
        cello_mappings_dp = join(self._data_dir, 'cellosaurus')
        
        makedirs(sources_dp) if not isdir(sources_dp) else None
        makedirs(cello_mappings_dp) if not isdir(cello_mappings_dp) else None

        cello_fp = self._download_cellosaurus_files(sources_dp)
        self._map_cello_names(cello_fp, join(cello_mappings_dp, 'cellosaurus_names.txt'))

    def _download_hpa_files(self, source_dir):
        """ Download HPA file

        Parameters
        ----------
        sources_dir : str
            the path to save the files

        Returns
        -------
        str
            the path to the HPA file
        """
        hpa_file_url = 'https://www.proteinatlas.org/download/proteinatlas.xml.gz'
        hpa_filepath = join(source_dir, "./proteinatlas.xml.gz")
        download_file_md5_check(hpa_file_url, hpa_filepath)
        return hpa_filepath

    def _parse_hpa_antibodies(self, hpa_fp, hpa_acc_fp, acc_hpa_fp):
        """ Map hpa tissues to cellosaurus cell lines

        Parameters
        ----------
        hap_fp : str
            the path to the hpa file
        hpa_acc_fp : str
            the path to output mappings from hpa antibodies to uniprot
        acc_hpa_fp : str
            the path to output mappings from uniprot to hpa antibodies
        """
        hpa_to_uniprot = defaultdict(set)
        uniprot_to_hpa = defaultdict(set)
        with gzip.open(hpa_fp) as xmlfile:
            for _, elem in tqdm(ET.iterparse(xmlfile), 'Processing HPA'):
                if elem.tag == 'entry' and len(elem) > 2:
                    id_elem = elem.find('identifier')
                    if len(id_elem) >= 1:
                        uniprot_ids = set()

                        xrefs = id_elem.findall('xref')
                        for xref in xrefs:
                            if xref.get('db') == 'Uniprot/SWISSPROT':
                                uniprot_ids.add(xref.get('id'))
                        if len(uniprot_ids) > 0:
                            antibodies = elem.findall('antibody')
                            if antibodies is not None:
                                for antibody in antibodies:
                                    ab_id = antibody.get('id')
                                    for uniprot_id in uniprot_ids:
                                        hpa_to_uniprot[ab_id].add(uniprot_id)
                                        uniprot_to_hpa[uniprot_id].add(ab_id)
                    elem.clear()


        with open(acc_hpa_fp, 'w', encoding='utf-8') as writer:
            for acc, antibodies in uniprot_to_hpa.items():
                if len(antibodies) > 0:
                    writer.write(f'{acc}\t{";".join(antibodies)}\n')

        with open(hpa_acc_fp, 'w', encoding='utf-8') as writer:
            for antibody, accs in hpa_to_uniprot.items():
                if len(accs) > 0:
                    writer.write(f'{antibody}\t{";".join(accs)}\n')

    def generate_hpa_mappings(self):
        """ Generate mappings for hpa"""
        sources_dp = self._source_dir
        
        hpa_mappings_dp = join(self._data_dir, 'hpa')
        uniprot_mappings_dp = join(self._data_dir, 'uniprot')
        makedirs(sources_dp) if not isdir(sources_dp) else None
        makedirs(hpa_mappings_dp) if not isdir(hpa_mappings_dp) else None
        makedirs(uniprot_mappings_dp) if not isdir(uniprot_mappings_dp) else None

        hpa_fp = self._download_hpa_files(sources_dp)
        self._parse_hpa_antibodies(
            hpa_fp,
            join(hpa_mappings_dp, 'hpa_to_acc.txt'),
            join(uniprot_mappings_dp, 'acc_to_hpa.txt')
        )

    def _clear_sources(self):
        """Delete the source files used to generate the mappings"""
        sources_dp = join(self._data_dir, 'sources')
        if exists(sources_dp):
            shutil.rmtree(sources_dp)
        
    def compress_mappings(self):
        filecount = 0
        mapping_sources = get_all_mappings_sources()
        for database, mappings in mapping_sources.items():
            filecount += len(mappings)

        t = tqdm(total=filecount, desc = 'Compressing mapping files')
        
        for root, dirs, files in walk(self._data_dir):
            for fp in files:
                src_fp = join(root, fp)
                dst_fp = join(root, fp+'.gz')
                with open(src_fp, 'rb') as f_in, gzip.open(dst_fp, 'wb') as f_out:
                    f_out.writelines(f_in)
                remove(src_fp)
                t.update(n=1)
        t.close()

    def verify_mappings(self):
        all_exist = True
        biodblinker_data = get_data_directory()
        mapping_sources = get_all_mappings_sources()
        for database, mappings in mapping_sources.items():
            for map_name, short_path in mappings.items():
                map_path = join(biodblinker_data, short_path)
                valid_map = exists(map_path) & isfile(map_path)
                if not exists(map_path):
                    print(f'{database}:{map_name} mapping file {map_path} does not exist')
                    all_exist = False
                elif not isfile(map_path):
                    print(f'{map_name} mapping file {map_path} is not a file')
                    all_exist = False
        
        if all_exist:
            print(f'All mappings file found')
        else:
            print(f'Could not find all mapping files')

    def generate_mappings(self, drugbank_user, drugbank_password, delete_sources=True):
        """ Generate mappings required by the biodblinker linkers

        Paramaters
        ----------
        drugbank_user : str
            drugbank username
        drugbank_password : str
            drugbank password
        """
        self.generate_cellosaurus_mappings()
        self.generate_hpa_mappings()
        self.generate_sider_mappings()
        self.generate_kegg_mappings()
        self.generate_drugbank_mappings(drugbank_user, drugbank_password)
        self.generate_uniprot_mappings()
        if delete_sources:
            self._clear_sources()
        self.compress_mappings()
