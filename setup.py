#!/usr/bin/env python3

import sys, os
from setuptools import setup

# ensure python version 3.5 or greater is used
if (sys.version_info.major + .1 * sys.version_info.minor) < 3.8:
	print('ERROR: please execute setup.py with python version >=3.8')
	sys.exit(1)

# get version string from seperate file
# https://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
# https://stackoverflow.com/questions/436198/what-is-an-alternative-to-execfile-in-python-3/437857#437857
VERSIONFILE="cov2seq/version.py"
with open(VERSIONFILE) as f:
	code = compile(f.read(), VERSIONFILE, 'exec')
	exec(code)
if not __version__:
	print("ERROR: unable to read version string from file {}".format(VERSIONFILE))
	exit()

DESCR = '''Comprehensive reports for Nanopore/Illumina/Sanger sequencing experiments of SARS CoV2 samples'''

# load long description from Markdown file
with open('README.md', 'rb') as readme:
	LONG_DESCR = readme.read().decode()

# check if defaults for user, host and dest are set for file transfer
setup_dir = os.path.dirname(os.path.abspath(__file__))

setup(name='cov2seq',
	  version=__version__,
	  description=DESCR,
	  long_description=LONG_DESCR,
	  url='http://github.com/MarkusHaak/cov2seq',
	  author='Markus Haak',
	  author_email='markus.haak@posteo.net',
	  license='GPL',
	  packages=['cov2seq'],
	  install_requires=['Biopython', 'vcf', 'numpy', 'pandas', 'matplotlib'],
	  include_package_data=True,
	  zip_safe=False,
	  entry_points={"console_scripts": ['cov2seq-update = cov2seq.main:update',
	  									'cov2seq-report = cov2seq.main:report']})
#	  scripts=['bin/cov2seq'])