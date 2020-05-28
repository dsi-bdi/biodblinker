# biodblinker
A library for linking entities of biological knowledge bases.

This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License.

# Installing
`pip install biodblinker`

# Installing from source

`python setup.py`

# Downloading mappings

When a biodblinker is initialized it verifies that all necessary mapping files are present and if not downloads the precompiled mappings

# Building the mapping files

It is also possible to generate the mappings from their sources

* Note this process will take several hours and requires a large ammount of disk space due to the size of the source files. The source files are removed once the mappings are generated

```
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

