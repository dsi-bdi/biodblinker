from biolink import KEGGLinker, GeneNameLinker, DrugBankLinker


kegg_db = KEGGLinker()
gene_name_db = GeneNameLinker()
drugbank_db = DrugBankLinker()

kegg_gene_ids = kegg_db.gene_ids
kegg_drug_ids = kegg_db.drug_ids

kegg_gene_ids_list = ["hsa:102723859", "hsa:64284", "hsa:80314", "hsa:79727", "hsa:106821730", "hsa:2909"]

uniprot_acc_list = kegg_db.convert_geneid_to_uniprot(kegg_gene_ids)
gene_names_list = kegg_db.convert_geneid_to_names(kegg_gene_ids)

single_names_list = [x[0] if len(x) > 0 else '' for x in gene_names_list]
gene_uniprots = gene_name_db.convert_gene_name_to_uniprot(single_names_list)
single_uniprots = [x[0] if len(x) > 0 else '' for x in gene_uniprots]
gene_uniprots = [d.split(";")[0] for d in single_uniprots]

print(f"= kegg gene ids: {kegg_gene_ids}")
print(f"= linked gene names: {gene_names_list}")

print(f"= [KG] uniprot ids: {uniprot_acc_list}")

print(f"= [GN] uniprot ids: {gene_uniprots}")
