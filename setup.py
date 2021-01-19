#!/usr/bin/env python3

import sys, os
from setuptools import setup
from glob import glob

# https://stackoverflow.com/questions/36187264
def binaries_directory():
    """Return the installation directory, or None"""
    if '--user' in sys.argv:
        paths = (site.getusersitepackages(),)
    else:
        py_version = '%s.%s' % (sys.version_info[0], sys.version_info[1])
        paths = (s % (py_version) for s in (
            sys.prefix + '/lib/python%s/dist-packages/',
            sys.prefix + '/lib/python%s/site-packages/',
            sys.prefix + '/local/lib/python%s/dist-packages/',
            sys.prefix + '/local/lib/python%s/site-packages/',
            '/Library/Python/%s/site-packages/',
        ))

    for path in paths:
        if os.path.exists(path):
            return path
    print('no installation path found', file=sys.stderr)
    return None


def submodule_data_files(submodule):
    data_files = []
    for dirpath, dirnames, filenames in os.walk(submodule):
        if filenames:
            install_dir = dirpath
            list_entry = (install_dir, [os.path.join(dirpath, fn) for fn in filenames])
            data_files.append(list_entry)
    return data_files

# ensure python version 3.7 or greater is used
if (sys.version_info.major + .1 * sys.version_info.minor) < 3.7:
    print('ERROR: please execute setup.py with python version >=3.7')
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

data_files = [('cov2seq', ['cfg/cov2seq.cfg'])]
data_files.extend(submodule_data_files('.git'))
data_files.extend(submodule_data_files('ncov'))

for d, files in data_files:
    print(d, files)

setup(name='cov2seq',
      version=__version__,
      description=DESCR,
      long_description=LONG_DESCR,
      url='http://github.com/MarkusHaak/cov2seq',
      author='Markus Haak',
      author_email='markus.haak@posteo.net',
      license='GPL',
      packages=['cov2seq'],
      install_requires=['Biopython', 'pyvcf', 'numpy', 'pandas', 'matplotlib', 'nextstrain-augur', 'Jinja2'],
      include_package_data=True,
      zip_safe=False,
      entry_points={"console_scripts": ['cov2seq-update = cov2seq.main:update',
                                        'cov2seq-report = cov2seq.main:report']},
      data_files=data_files,
      package_data={'cov2seq': ['templates/report.html']})

site_packages_dir = binaries_directory()
if site_packages_dir:
    pkg_dir = "cov2seq-{}-py{}.egg".format(__version__, sys.version_info.major + .1 * sys.version_info.minor)
    default_cfg = os.path.join(site_packages_dir, pkg_dir, "cov2seq", "cov2seq.cfg")
    user_cfg = os.path.join(os.path.expanduser('~'), "cov2seq.cfg")
    print(default_cfg)
    if os.path.exists(default_cfg):
        print("\nINFO: default configuration file is found under {}".format(default_cfg))
        print("If (partial) configurations are given in user config file {}, they are prioritized.".format(user_cfg))
    else:
        print("\nWARNING: default configuration file not found under {}".format(default_cfg))
        print("Please supply all configurations found in configuration file \"{}\" in a user configuration file under {}".format("cov2seq.cfg", user_cfg))

