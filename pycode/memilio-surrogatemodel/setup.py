import os

from setuptools import find_packages, setup

__version__ = '0.1.0'

setup(
    name='memilio-surrogatemodel', version=__version__, author='DLR-SC',
    author_email='henrik.zunker@dlr.de',
    maintainer_email='martin.kuehn@dlr.de',
    url='https://github.com/DLR-SC/memilio',
    description='Part of MEmilio project, implementation of surrogate models for the existing models in MEmilio.',
    packages=find_packages(
        where=os.path.dirname(os.path.abspath(__file__))),
    install_requires=[
        # smaller pandas versions contain a bug that sometimes prevents reading
        # some excel files (e.g. population or twitter data)'numpy',
        'pandas>=1.2.2',
        'progress',
        'numpy>=1.22',  # smaller numpy versions cause a security issue
        'tensorflow',
        'matplotlib', ],
    extras_require={'dev': [
        # smaller pyfakefs versions use deprecated functions for matplotlib versions >=3.4
        'pyfakefs>=4.2.1',
        'coverage',
    ], },
    long_description='', test_suite='memilio.surrogatemodel_test',)
