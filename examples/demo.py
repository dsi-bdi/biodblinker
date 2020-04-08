from biolink import (KEGGLinker, GeneNameLinker, DrugBankLinker,
                     UniprotLinker, HPALinker, CellosaurusLinker)


kegg_db = KEGGLinker()
gene_name_db = GeneNameLinker()
drugbank_db = DrugBankLinker()
uniprot_db = UniprotLinker()
hpa_db = HPALinker()
cl_db = CellosaurusLinker()

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

print()
cells = ['CVCL_0028', 'CVCL_0025', 'CVCL_Y019']
tissues = ['Lung', 'Brain']
cell_tissues = cl_db.convert_cell_line_to_hpa(cells)
cell_names = cl_db.convert_cell_line_to_name(cells)
tissue_cells = hpa_db.convert_hpa_to_cellosaurus(tissues)

print()
print('Cellosaurus cell line ids mapped to hpa tissues')
for index, cell in enumerate(cells):
    print(f'{cell} linked to {cell_tissues[index]}')

print()
print('Cellosaurus cell line ids mapped to names')
for index, cell in enumerate(cells):
    print(f'{cell} linked to {cell_names[index]}')



print()
print('HPA tissues mapped to Cellosaurus cell line ids')
for index, tissue in enumerate(tissues):
    print(f'{tissue} linked to {tissue_cells[index]}')
