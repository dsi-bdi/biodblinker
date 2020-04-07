from os import makedirs
import gzip
import os.path
from biolink.fileio import *
from biolink.config import file_open
from tqdm import tqdm
import requests
from collections import defaultdict
from zipfile import ZipFile
import xml.etree.ElementTree as ET
# ================================================================================
# Loading UNIPROT mapping for human accession ids
# ================================================================================
target_databases = {"Gene_Name", "HPA", "PharmGKB", "STRING", "HGNC", "BioGrid", "RefSeq", "UniParc", "KEGG", "Ensembl",
                    "MIM", "MINT", "UniProtKB-ID", "EMBL-CDS", "DIP", "ProteomicsDB", "UniRef100", "GI", "GeneDB",
                    "EMBL", "ChEMBL", "KO", "MIM", "Gene_Synonym"}

mapping_data_root_dp = "../data"
uniprot_mappings_dp = join(mapping_data_root_dp, "uniprot")
mapping_data_srcs_dp = join(mapping_data_root_dp, "sources")
makedirs(mapping_data_srcs_dp) if not isdir(mapping_data_srcs_dp) else None
makedirs(uniprot_mappings_dp) if not isdir(uniprot_mappings_dp) else None

uniprot_mapping_file_url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz"
uniprot_mapping_filepath = join(mapping_data_srcs_dp, "./uniprot_hsa_mappings.dat.gz")
download_file_md5_check(uniprot_mapping_file_url, uniprot_mapping_filepath)

uniprot_map_fd = file_open(uniprot_mapping_filepath)
db_mapping_dictionaries_dict = {dbname: dict() for dbname in target_databases}

for line in tqdm(uniprot_map_fd, desc="Processing UNIPROT mapping file."):
    uniprot_acc, db_tag, db_id = line.strip().split("\t")
    # change "Gene_Synonym" to "Gene_Name" to treat both equally

    if db_tag in target_databases:
        if uniprot_acc not in db_mapping_dictionaries_dict[db_tag]:
            db_mapping_dictionaries_dict[db_tag][uniprot_acc] = [db_id]
        else:
            db_mapping_dictionaries_dict[db_tag][uniprot_acc].append(db_id)

for db_name in tqdm(db_mapping_dictionaries_dict, desc="Exporting uniprot target databases mappings"):
    db_map_fp = join(uniprot_mappings_dp, f"acc_to_{db_name.lower()}.txt")
    db_map_fd = open(db_map_fp, "w")
    for key, val in db_mapping_dictionaries_dict[db_name].items():
        val_txt = ";".join(val)
        db_map_fd.write(f"{key}\t{val_txt}\n")
    db_map_fd.close()

# ================================================================================
# KEGG
# ================================================================================
class SetWriter:
    """
    Utility class for writing DrugBank statements
    Enforces uniqueness of statements between written between flushes
    Set clear_on_flush to false to enforce uniquness for on all writes
    (should not be set for very large datasets)
    """
    def __init__(self, path):
        """
        Initialize a new SetWriter

        Parameters
        ----------
        """
        self._lines = []
        self._lineset = set()
        if not os.path.exists(os.path.dirname(path)):
            makedirs(os.path.dirname(path))
        self._fd = open(path, 'w')
        self._clear_on_flush = True
        self._closed = False

    @property
    def clear_on_flush(self):
        return self._clear_on_flush

    @clear_on_flush.setter
    def clear_on_flush(self, value):
        self._clear_on_flush = value

    def write(self, line):
        if self._closed:
            raise ValueError('I/O operation on closed file')

        if line in self._lineset:
            return
        self._lineset.add(line)
        self._lines.append(line)

    def flush(self):
        if self._closed:
            raise ValueError('I/O operation on closed file')
        self._fd.writelines(self._lines)
        self._lines = []
        if self._clear_on_flush:
            self._lineset = set()

    def close(self):
        if len(self._lines) > 0:
            self.flush()
        self._lineset = set()
        self._fd.close()


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
    'compound':
    [
        '3dmet', 'hsdb', 'hmdb', 'nikkaji', 'chembl',
        'knapsack', 'pubchem', 'chebi', 'pdb-ccd',
        'lipidbank', 'lipidmaps', 'massbank'
    ],
    'hsa':
    [
        'uniprot', 'ensembl'
    ],
    'pathway': [],
    'brite': [],
    'module': [],
    'ko': [],
    'glycan': [],
    'reaction': [],
    'enzyme': [],
    'network': [],
    'variant': [],
    'dgroup': [],
    'environ': []
}


def map_kegg_names(database):
    id_set = set()
    database_name = database
    if database == 'hsa':
        database_name = 'gene'
    writer = SetWriter(f'../data/kegg/{database_name}_names.txt')
    resp = requests.get(f'http://rest.kegg.jp/list/{database}')

    if resp.ok:
        for line in resp.iter_lines(decode_unicode=True):
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                (source_str, name_str) = parts[0:2]
                pre, did = source_str.split(':')
                if did[0].isdigit():
                    did = source_str
                id_set.add(did)

                names = set()
                for name in name_str.split(';'):
                    name = name.strip()
                    if database == 'drug' or database == 'environw':
                        index = name.rfind('(')
                        if index > -1:
                            name = name[:index].strip()

                    names.add(name)
                writer.write(f'{did}\t{";".join(names)}\n')
            else:
                source_str = parts[0]
                pre, did = source_str.split(':')
                if did[0].isdigit():
                    did = source_str
                id_set.add(did)
                writer.write(f'{did}\t-\n')

    writer.close()
    return id_set


def link_kegg_db(database, target):
    mapping_dict = defaultdict(set)
    resp = requests.get(f'http://rest.genome.jp/link/{target}/{database}')
    is_empty = True
    if resp.ok:
        for line in resp.iter_lines(decode_unicode=True):
            is_empty = False
            source_str, target_str, _ = line.strip().split('\t')
            pre, did = source_str.split(':')
            if did[0].isdigit():
                did = source_str

            _, tid = target_str.split(':')
            if target == 'sider':
                sider_base = 'CID100000000'
                tid = sider_base[:len(sider_base)-len(tid)]+tid
            mapping_dict[did].add(tid)

    if is_empty:
        print(f'Unable to link {database} to {target}')
    return mapping_dict


for source, targets in tqdm(kegg_source_targets.items(), 'Processing KEGG mappings'):
    source_ids = map_kegg_names(source)
    for target in targets:
        mapping_dict = link_kegg_db(source, target)
        if len(mapping_dict) > 0:
            source_name = source
            if source == 'hsa':
                source_name = 'gene'

            reverse_map = defaultdict(set)
            writer = SetWriter(f'../data/kegg/{source_name}_{target}.txt')
            for sid in source_ids:
                map_ids = mapping_dict[sid]
                if len(map_ids) == 0:
                    map_ids.add('-')
                else:
                    for mapped_id in map_ids:
                        reverse_map[mapped_id].add(sid)
                writer.write(f'{sid}\t{";".join(map_ids)}\n')
            writer.close()

            writer = SetWriter(f'../data/{target}/{target}_kegg_{source_name}.txt')
            for tid, kids in reverse_map.items():
                writer.write(f'{tid}\t{";".join(kids)}\n')
            writer.close()

# ==============================================================================
# DRUGBANK
# ==============================================================================
drugbank_targets = {
    'BindingDB': 'bindingdb',
    'ChEBI': 'chebi',
    'ChEMBL': 'chembl',
    'ChemSpider': 'chemspider',
    'Drugs Product Database (DPD)': 'dpd',
    'IUPHAR': 'iuphar',
    'KEGG Compound': 'kegg_compound',
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

def map_drugbank():
    id_set = set()
    drugbank_targets_map = {}
    drugbank_targets_map['names'] = defaultdict(set)
    drugbank_targets_map['ids'] = defaultdict(set)
    for res, target in drugbank_targets.items():
        drugbank_targets_map[target] = defaultdict(set)

    ns = {'db':'http://www.drugbank.ca'}
    with ZipFile('../data/sources/drugbank_all_full_database.xml.zip', 'r') as dbzip:
        with dbzip.open('full database.xml', force_zip64=True) as xmlfile:
            for event, elem in tqdm(ET.iterparse(xmlfile), 'Processing Drugbank mapping'):
                # Check the length of the drug element as pathways also contain drug elements
                if elem.tag =='{http://www.drugbank.ca}drug' and len(elem) > 2:
                    drug_id_elem = elem.find('./db:drugbank-id[@primary="true"]', ns)
                    did = drug_id_elem.text
                    
                    id_set.add(did)
                    name_elem = elem.find('./db:name', ns)
                    names = set()
                    if name_elem is not None:
                        drugbank_targets_map['names'][did].add(name_elem.text.strip())
                        
                    synonyms = elem.find('./db:synonyms', ns)
                    if synonyms is not None:
                        for synonym in synonyms:
                            drugbank_targets_map['names'][did].add(synonym.text.strip())

                    
                    id_elems = elem.findall('./db:drugbank-id', ns)
                    ids = set()
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
                        else:
                            print(source.text)

                    elem.clear()
    return id_set, drugbank_targets_map


drugbank_ids, drugbank_targets_map = map_drugbank()
for target, target_map in tqdm(drugbank_targets_map.items(), 'Exporting Drugbank mappings'):
    if len(target_map) == 0:
        continue
    reverse_map = defaultdict(set)
    writer = SetWriter(f'../data/drugbank/drugbank_{target}.txt')
    for did in sorted(drugbank_ids):
        mapped = target_map[did]
        if len(mapped) == 0:
            mapped.add('-')
        else:
            for mapped_id in mapped:
                reverse_map[mapped_id].add(did)
        writer.write(f'{did}\t{";".join(mapped)}\n')
    writer.close()

    writer = SetWriter(f'../data/{target}/{target}_drugbank.txt')
    for tid, dids in reverse_map.items():
        writer.write(f'{tid}\t{";".join(dids)}\n')
    writer.close()
    
# ==============================================================================
# SIDER
# ==============================================================================
sider_names_file = 'http://sideeffects.embl.de/media/download/drug_names.tsv'
sider_names_filepath = join(mapping_data_srcs_dp, "./sider_drug_names.tsv")
download_file_md5_check(sider_names_file, sider_names_filepath)
stitch_mapping_file = 'http://stitch.embl.de/download/chemical.sources.v5.0.tsv.gz'
stitch_mapping_filepath = join(mapping_data_srcs_dp, "./stitch_mapping.tsv.gz")
download_file_md5_check(stitch_mapping_file, stitch_mapping_filepath)

sider_targets = {
    'BindingDB': 'bindingdb',
    'ChEBI': 'chebi',
    'ChEMBL': 'chembl',
    'DrugBank': 'drugbank',
    'KEGG': 'kegg',
    'pc': 'pubchem_compound',
    'ps': 'pubchem_substance',
    'ATC': 'atc'
}


def map_sider():
    sider_ids = set()
    sider_target_map = {}
    sider_target_map['names'] = defaultdict(set)
    for res, target in sider_targets.items():
        sider_target_map[target] = defaultdict(set)

    with open('../data/sources/sider_drug_names.tsv', 'r') as names_fd:
        for line in names_fd:
            sid, name = line.strip().split('\t')
            sider_target_map['names'][sid].add(name)
            sider_ids.add(sid)

    with gzip.open('../data/sources/stitch_mapping.tsv.gz', 'rt') as stitch_fd:
        next(stitch_fd)
        for line in tqdm(stitch_fd, 'Processing Sider mapping'):
            if line.startswith('#'):
                continue
            (mole, subs, target, target_id) = line.strip().split('\t')
            sid = 'CID1'+mole[4:]
            if target in sider_targets:
                target_name = sider_targets[target]
                sider_target_map[target_name][sid].add(target_id)

    return sider_ids, sider_target_map


sider_ids, sider_target_map = map_sider()
for target, target_map in tqdm(sider_target_map.items(), 'Exporting Sider mapping'):
    if len(target_map) == 0:
        continue

    reverse_map = defaultdict(set)
    writer = SetWriter(f'../data/sider/sider_{target}.txt')
    for sid in sorted(sider_ids):
        mapped = target_map[sid]
        if len(mapped) == 0:
            mapped.add('-')
        else:
            for mapped_id in mapped:
                reverse_map[mapped_id].add(sid)
        writer.write(f'{sid}\t{";".join(mapped)}\n')
    writer.close()

    writer = SetWriter(f'../data/{target}/{target}_sider.txt')
    for tid, sids in reverse_map.items():
        writer.write(f'{tid}\t{";".join(sids)}\n')
    writer.close()

