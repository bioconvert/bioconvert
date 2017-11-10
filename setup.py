# -*- coding: utf-8 -*-
import os
from setuptools import setup, find_packages

_MAJOR = 0
_MINOR = 0
_MICRO = 6
version = '%d.%d.%d' % (_MAJOR, _MINOR, _MICRO)
release = '%d.%d' % (_MAJOR, _MINOR)

metainfo = {
    'authors': {
        'Cokelaer': ('Thomas Cokelaer', 'thomas.cokelaer@pasteur.fr'),
        },
    'version': version,
    'license': 'BSD',
    'download_url': ['http://pypi.python.org/pypi/bioconvert'],
    'url': ['http://pypi.python.org/pypi/bioconvert'],
    'description': 'convert between bioinformatics formats',
    'platforms': ['Linux', 'Unix', 'MacOsX', 'Windows'],
    "keywords": ["NGS", "bam2bed", "fastq2fasta", "bam2sam"],
    'classifiers': [
          'Development Status :: 1 - Planning',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Topic :: Software Development :: Libraries :: Python Modules',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Information Analysis',]
    }


with open('README.rst') as f:
    readme = f.read()

requirements = open("requirements.txt").read().split()

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd:
    # mock, pillow, sphinx, sphinx_rtd_theme installed on RTD
    # but we also need numpydoc and sphinx_gallery
    extra_packages = ["numpydoc", "sphinx_gallery"]
    requirements += extra_packages

# pygatb name on pypi is pygatb (lower case) but is not in bioconda.
# Since the requirements.txt is used both by this file and .travis (bioconda),
# travis will fail. So, we must hard-code remove pygatb from the
# requirements.txt Consequently, we have to hard-code it in this script.
requirements += ["pygatb"]

setup(
    name='bioconvert',
    version=version,
    maintainer=metainfo['authors']['Cokelaer'][0],
    maintainer_email=metainfo['authors']['Cokelaer'][1],
    author='The bioconvert Contributors',
    author_email=metainfo['authors']['Cokelaer'][1],
    long_description=readme,
    keywords=metainfo['keywords'],
    description=metainfo['description'],
    license=metainfo['license'],
    platforms=metainfo['platforms'],
    url=metainfo['url'],
    download_url=metainfo['download_url'],
    classifiers=metainfo['classifiers'],
    zip_safe=False,
    packages=find_packages(),
    install_requires=requirements,

    # This is recursive include of data files
    exclude_package_data={"": ["__pycache__"]},

    package_data={
        '': ['*.csv'],
        'bioconvert.data': ['*'],
        },

    entry_points={
        'console_scripts': [
           'bioconvert=bioconvert.scripts.converter:main',
           'bioconvert_init=bioconvert.scripts.init_convert:main'
        ]
    }

    )
