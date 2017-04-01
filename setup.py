"""Software pipeline for end-to-end analysis of dual-CRISPR screens

This software implements a Jupyter Notebook-based pipeline for
automated end-to-end analysis of cell-based dual-CRISPR screens.
It has been developed for the lab of Dr. Prashant Mali by Amanda
Birmingham and Roman Sasik of  the Center for Computational Biology
and Bioinformatics at the University of California at San Diego.
"""

# Much of the content of this file is copied from the
# setup.py of the (open-source) PyPA sample project at
# https://github.com/pypa/sampleproject/blob/master/setup.py

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='dual_crispr',

    # Versions should comply with PEP440.
    version='0.3.0',

    description='Software pipeline for end-to-end analysis of dual-CRISPR screens',
    long_description=long_description,

    # The project's main homepage.
    url="https://github.com/ucsd-ccbb/mali-dual-crispr-pipeline",

    # Author details
    author='The Center for Computational Biology and Bioinformatics',
    author_email='abirmingham@ucsd.edu',

    # Choose your license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.'

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language:: Python:: 3:: Only',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    # What does your project relate to?
    keywords='crispr high-throughput screening bioinformatics',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed, although this can be overridden with a requirements.txt file
    install_requires=['ccbb_pyutils', 'cutadapt',
                      'matplotlib', 'numpy', 'pandas', 'rpy2'],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    extras_require={
        'dev': ['check-manifest'],
        'test': ['coverage'],
    },

    package_data={
        'dual_crispr': ['distributed_files/notebooks/*.ipynb',
                        'distributed_files/library_definitions/*.txt',
                        'distributed_files/test_data/TestRun1/*.fastq',
                        'distributed_files/test_data/TestRun2/*.txt',
                        'distributed_files/test_data/TestRun3/*.txt',
                        'distributed_files/test_data/TestRun3withError/*.txt']
    },

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'set_up_dual_crispr=dual_crispr.run_set_up_pipeline:main',
            'count_dual_crispr=dual_crispr.run_counting:main',
            'score_dual_crispr=dual_crispr.run_scoring:main'
        ]
    }
)