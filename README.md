# biodblinker
A library for linking entities of biological knowledge bases.

This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License.

# Installing
`pip install biodblinker`

# Installing from source

`python setup.py`

# Usage
```python
from biodblinker import UniprotLinker

uniprot_linker = UniprotLinker()

# Get list of all included uniprot accessions
uniprot_accessions = uniprot_linker.uniprot_ids

select_accessions = ['P31946', 'P62258', 'Q04917']

# Get the list of names for each accession in select_accessions
select_names = uniprot_linker.convert_uniprot_to_names(select_accessions)

# Get the list of kegg gene ids for each accession in select_accessions
select_genes = uniprot_linker.convert_uniprot_to_kegg(select_accessions)

```

# Use Case - Linking uniprot proteins and mesh diseases via KEGG

```python
import requests
from biodblinker import KEGGLinker

linker = KEGGLinker()
unique_pairs = set()

url = 'http://rest.kegg.jp/link/hsa/disease'
resp = requests.get(url)

if resp.ok:
    for line in resp.iter_lines(decode_unicode=True):
        kegg_disease, kegg_gene = line.strip().split('\t')
        # strip the prefix from the disease
        kegg_disease = kegg_disease.split(':')[1]
        
        # prefix is retained for genes as the ids an numeric
        uniprot_protein = linker.convert_geneid_to_uniprot([kegg_gene])
        mesh_disease = linker.convert_disease_to_mesh([kegg_disease])
        if len(uniprot_protein[0]) == 0:
            continue
        if len(mesh_disease[0]) == 0:
            continue
        for protein in uniprot_protein[0]:
            for disease in mesh_disease[0]:
                unique_pairs.add((protein, disease))

for protein, disease in unique_pairs:
    print(f'{protein}\tRELATED_DISEASE\t{disease}')

```

# Downloading mappings

When a biodblinker is initialized it verifies that all necessary mapping files are present and if not downloads the precompiled mappings

# Building the mapping files

It is also possible to generate the mappings from their sources

* Note this process will take several hours and requires a large ammount of disk space due to the size of the source files. The source files are removed once the mappings are generated

```python
import biodblinker

gen = biodblinker.MappingGenerator()
gen.generate_mappings(<drugbank_username>, <drugbank_password>)
```


# Mapping sources and licenses
BioDBLinker uses multiple sources to generate the mappings. BioDBLinker must be used in compliance with these licenses and citation policies where applicable

| Source Database                                    | License Type | URL                                                |
|----------------------------------------------------|--------------|----------------------------------------------------|
| [UniProt](https://www.uniprot.org)                 | CC BY 4.0    | https://www.uniprot.org/help/license               |
| [Drugbank](https://www.drugbank.ca/)               | CC BY NC 4.0 | https://www.drugbank.ca/legal/terms_of_use         |
| [KEGG](https://www.genome.jp/kegg/)                | Custom       | https://www.kegg.jp/kegg/legal.html                |
| [Sider](http://sideeffects.embl.de/)               | CC BY-NC-SA  | http://sideeffects.embl.de/about/                  |
| [Stitch](http://stitch.embl.de/)                   | CC BY 4.0    | http://stitch.embl.de/cgi/download.pl              |
| [HPA](https://www.proteinatlas.org/)               | CC BY SA 3.0 | https://www.proteinatlas.org/about/licence         |
| [Cellosaurus](https://web.expasy.org/cellosaurus/) | CC BY 4.0    | https://web.expasy.org/cgi-bin/cellosaurus/faq#Q22 |

