# biolink
A library for linking entities of biological knowledge bases.

# building
`python setup.py sdist`

`conda build ./conda_recipe`

`conda install --use-local biolink`

# package for cross platform
After building for conda

`conda convert -p all <path to built package>`

\<path to build package\> is output in the logs of conda build

# building the mapping files

```
import biolink

gen = biolink.MappingGenerator()
gen.generate_mappings(<drugbank_username>, <drugbank_password>)
```