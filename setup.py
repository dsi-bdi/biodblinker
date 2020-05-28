from setuptools import setup
from biodblinker import __version__, __name__, __description__

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
   name=__name__,
   version=__version__,
   description=__description__,
   author='INSIGHT Centre for Data Analytics',
   author_email='sameh.kamal@insight-centre.org',
   long_description=long_description,
   long_description_content_type='text/markdown',
   url='https://github.com/dsi-bdi/biodblinker',
   packages=[__name__],
   include_package_data=True,
   install_requires=[
      'tqdm',
      'requests'
   ],
   
)
