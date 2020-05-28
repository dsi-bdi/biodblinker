from os import makedirs
import gzip
from biodblinker.fileio import *
from biodblinker.config import file_open
from tqdm import tqdm

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

uniprot_mapping_file_url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz"
uniprot_mapping_filepath = join(mapping_data_srcs_dp, "./uniprot_hsa_mappings.dat.gz")
download_file_md5_check(uniprot_mapping_file_url, uniprot_mapping_filepath)

uniprot_map_fd = file_open(uniprot_mapping_filepath)

db_mapping_dictionaries_dict = {dbname: dict() for dbname in target_databases}
db_mapping_dictionaries_reverse_dict = {dbname: dict() for dbname in target_databases}

for line in tqdm(uniprot_map_fd, desc="Processing UNIPROT mapping file."):
    uniprot_acc, db_tag, db_id = line.strip().split("\t")
    # change "Gene_Synonym" to "Gene_Name" to treat both equally

    if db_tag in target_databases:
        if uniprot_acc not in db_mapping_dictionaries_dict[db_tag]:
            db_mapping_dictionaries_dict[db_tag][uniprot_acc] = [db_id]
        else:
            db_mapping_dictionaries_dict[db_tag][uniprot_acc].append(db_id)

        if db_id not in db_mapping_dictionaries_reverse_dict[db_tag]:
            db_mapping_dictionaries_reverse_dict[db_tag][db_id] = [uniprot_acc]
        else:
            db_mapping_dictionaries_reverse_dict[db_tag][db_id].append(uniprot_acc)

# ================================================================================
# Export mappings from UNIPROT to other databases
# ================================================================================
for db_name in tqdm(db_mapping_dictionaries_dict, desc="Exporting uniprot target databases mappings"):
    db_map_fp = join(uniprot_mappings_dp, f"acc_to_{db_name.lower()}.txt")
    db_map_fd = open(db_map_fp, "w")
    for key, val in db_mapping_dictionaries_dict[db_name].items():
        val_txt = ";".join(val)
        db_map_fd.write(f"{key}\t{val_txt}\n")
    db_map_fd.close()

# ================================================================================
# Export mappings from databases to UNIPROT
# ================================================================================

# add gene synonyms to gene name database entry in the reverse maps dictionary and remove the synonyms entry
db_mapping_dictionaries_reverse_dict["Gene_Name"].update(db_mapping_dictionaries_reverse_dict["Gene_Synonym"].copy())
del db_mapping_dictionaries_reverse_dict["Gene_Synonym"]

# export to files in respective dictionaries
for db_name in tqdm(db_mapping_dictionaries_reverse_dict, desc="Exporting uniprot target databases reverse mappings"):
    db_mappings_dp = join(mapping_data_root_dp, db_name.lower())
    makedirs(db_mappings_dp) if not isdir(db_mappings_dp) else None
    db_rev_map_fp = join(db_mappings_dp, f"{db_name.lower()}_to_uniprot.txt")
    db_rev_map_fd = open(db_rev_map_fp, "w")
    for key, val in db_mapping_dictionaries_reverse_dict[db_name].items():
        val_txt = ";".join(val)
        db_rev_map_fd.write(f"{key}\t{val_txt}\n")
    db_rev_map_fd.close()
# ================================================================================




