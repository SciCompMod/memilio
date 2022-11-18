import os
from setuptools import setup, find_packages

__version__ = '0.1.0'

setup(
    name='memilio-surrogatemodel', version=__version__, author='DLR-SC',
    author_email='henrik.zunker@dlr.de',
    maintainer_email='martin.kuehn@dlr.de',
    url='https://github.com/DLR-SC/memilio',
    description='Part of MEmilio project, implementation of surrogate models for the existing models in MEmilio.',
    packages=find_packages(
        where=os.path.dirname(os.path.abspath(__file__))),
    install_requires=[],
    extras_require={'dev': ['numpy >= 1.22'], },
    long_description='', test_suite='memilio.surrogatemodel_test',)
