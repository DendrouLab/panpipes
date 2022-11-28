# setup.py script adapted from COMBAT sc_pipelines repo, originally by Adam Cribbs

import sys
import os
import re
import setuptools
from setuptools import setup, find_packages, Extension

from packaging.version import Version
if Version(setuptools.__version__) < Version('1.1'):
    print("Version detected:", Version(setuptools.__version__))
    raise ImportError(
        "sc_pipelines requires setuptools 1.1 higher")

########################################################################
########################################################################
IS_OSX = sys.platform == 'darwin'

########################################################################
########################################################################
# collect version
# not sure if this is necessary
sys.path.insert(0, "panpipes")
import version
version = version.__version__

###############################################################
###############################################################
# Define dependencies
#
major, minor1, minor2, s, tmp = sys.version_info

if major < 3:
    raise SystemExit("""Requires Python 3 or later.""")

sc_pipelines_packages = find_packages()
# print(sc_pipelines_packages)
sc_pipelines_package_dirs = {'pipelines': 'pipelines'}

##########################################################
##########################################################
# Classifiers
classifiers = """
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

setup(
    # package information
    name='panpipes',
    version=version,

    description='panpipes : single-cell pipelines',
    author='Charlotte Rich-Griffin',
    author_email='charlotte.rich@well.ox.ac.uk',
    license="MIT",
    platforms=["any"],
    keywords="computational genomics",
    long_description='''panpipes: Dendrou group single-cell pipelines''',

    classifiers=[_f for _f in classifiers.split("\n") if _f],
    url="",
    # package contents
    packages=['panpipes', 'panpipes/funcs', 'python', 'R', 'resources'],
    # packages=sc_pipelines_packages,
    # package_dir=sc_pipelines_package_dirs,
    include_package_data=True,
    entry_points={
        "console_scripts": ["panpipes = panpipes.entry:main"]
    },
    # other options
    zip_safe=False,
    test_suite="tests",
)
