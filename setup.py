from setuptools import setup
from biolink import __version__, __name__, __description__

setup(
   name=__name__,
   version=__version__,
   description=__description__,
   author='INSIGHT Centre for Data Analytics',
   author_email='sameh.kamal@insight-centre.org',
   packages=[__name__],
   include_package_data=True,
   install_requires=[
      'tqdm',
      'requests'
                     ],
)
