# -*- coding: utf-8 -*-
import os
from setuptools import setup, find_packages

_MAJOR = 0
_MINOR = 6
_MICRO = 2
version = "%d.%d.%d" % (_MAJOR, _MINOR, _MICRO)
release = "%d.%d" % (_MAJOR, _MINOR)

metainfo = {
    "authors": {"Cokelaer": ("Thomas Cokelaer", "thomas.cokelaer@pasteur.fr"),},
    "version": version,
    "license": "GPLv3",
    "download_url": "https://github.com/bioconvert/bioconvert",
    "url": "https://github.com/bioconvert/bioconvert",
    "description": "convert between bioinformatics formats",
    "platforms": ["Linux", "Unix", "MacOsX", "Windows"],
    "keywords": ["NGS", "bam2bed", "fastq2fasta", "bam2sam"],
    "classifiers": [
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Information Analysis",
    ],
}


with open("README.rst") as f:
    readme = f.read()

requirements = open("requirements.txt").read().split()
on_rtd = os.environ.get("READTHEDOCS", None) == "True"
if on_rtd:
    # mock, pillow, sphinx, sphinx_rtd_theme installed on RTD
    # but we also need numpydoc and sphinx_gallery
    extra_packages = ["numpydoc", "sphinx_gallery"]
    requirements += extra_packages


setup(
    name="bioconvert",
    version=version,
    maintainer=metainfo["authors"]["Cokelaer"][0],
    maintainer_email=metainfo["authors"]["Cokelaer"][1],
    author="The bioconvert Contributors",
    author_email=metainfo["authors"]["Cokelaer"][1],
    long_description=readme,
    long_description_content_type="text/x-rst",
    keywords=metainfo["keywords"],
    description=metainfo["description"],
    license=metainfo["license"],
    platforms=metainfo["platforms"],
    url=metainfo["url"],
    download_url=metainfo["download_url"],
    classifiers=metainfo["classifiers"],
    zip_safe=False,
    packages=find_packages(),
    install_requires=requirements,
    extras_require={
        "testing": [
            "pytest",
            "pytest-cov",
            "pytest-env",
            "pytest-xdist",
            "pytest-mock",
            "pytest-timeout",
            "pytest-runner",
            "pytest-rerunfailures",
            "mock",
            "coveralls",
        ],
        "doc": [
            "pillow",
            "sphinx",
            "sphinx_rtd_theme",
            "sphinx_gallery",
            "numpydoc",
            "pygraphviz",
        ],
    },
    # This is recursive include of data files
    exclude_package_data={"": ["__pycache__"]},
    package_data={
        "": ["*.csv", "*.sh"],
        "bioconvert.data": ["*"],
        "bioconvert.misc": ["*"],
    },
    entry_points={
        "console_scripts": [
            "bioconvert=bioconvert.scripts.converter:main",
            "bioconvert_init=bioconvert.scripts.init_convert:main",
            "bioconvert_stats=bioconvert.scripts.stats:main",
        ]
    },
)
