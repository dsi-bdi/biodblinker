from biolink import KEGGLinker, GeneNameLinker, DrugBankLinker, UniprotLinker


kegg_db = KEGGLinker()
gene_name_db = GeneNameLinker()
drugbank_db = DrugBankLinker()
uniprot_db = UniprotLinker()

kegg_gene_ids = kegg_db.gene_ids
kegg_drug_ids = kegg_db.drug_ids

kegg_gene_ids_list = ["hsa:102723859", "hsa:64284", "hsa:80314", "hsa:79727", "hsa:106821730", "hsa:2909"]

uniprot_acc_list = kegg_db.convert_geneid_to_uniprot(kegg_gene_ids_list)
gene_names_list = kegg_db.convert_geneid_to_names(kegg_gene_ids_list)

single_names_list = [x[0] if len(x) > 0 else '' for x in gene_names_list]
gene_uniprots = gene_name_db.convert_gene_name_to_uniprot(single_names_list)
single_uniprots = [x[0] if len(x) > 0 else '' for x in gene_uniprots]
gene_uniprots = [d.split(";")[0] for d in single_uniprots]

uniprot_acc_ids_list = [acc for linked_accs in uniprot_acc_list for acc in linked_accs]
linked_kegg_gene_ids = uniprot_db.convert_uniprot_to_kegg(uniprot_acc_ids_list)
linked_gene_names = uniprot_db.convert_uniprot_to_gene_name(uniprot_acc_ids_list)

print(f"= kegg gene ids: {kegg_gene_ids}")

print(f"= linked gene names: {gene_names_list}")

print(f"= [KG] uniprot ids: {uniprot_acc_list}")

print(f"= [GN] uniprot ids: {gene_uniprots}")

print('Kegg genes mapped to uniprot accs')
for index, kegg_gene in enumerate(kegg_gene_ids_list):
    print(f'{kegg_gene} linked to {uniprot_acc_list[index]}')

print()
print('Uniprot accs mapped to Kegg genes')
for index, acc in enumerate(uniprot_acc_ids_list):
    print(f'{acc} linked to {linked_kegg_gene_ids[index]}')

print()
print('Uniprot accs mapped to gene names')
for index, acc in enumerate(uniprot_acc_ids_list):
    print(f'{acc} linked to {linked_gene_names[index]}')



