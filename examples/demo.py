from biolink import KEGGLinker


kegg_db = KEGGLinker()

kegg_gene_ids = kegg_db.gene_ids
kegg_drug_ids = kegg_db.drug_ids

kegg_gene_ids_list = ["hsa:102723859", "hsa:64284", "hsa:80314", "hsa:79727", "hsa:106821730", "hsa:2909"]

uniprot_acc_list = kegg_db.convert_geneid_to_uniprot(kegg_gene_ids)
gene_names_list = kegg_db.convert_geneid_to_names(kegg_gene_ids)

print(f"= kegg gene ids: {kegg_gene_ids}")
print(f"= linked uniprot ids: {uniprot_acc_list}")
print(f"= linked gene names: {gene_names_list}")
