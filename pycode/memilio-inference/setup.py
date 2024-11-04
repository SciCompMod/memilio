import os

from setuptools import find_packages, setup

__version__ = '1.0.0'

setup(
    name='memilio-inference', version=__version__, author='DLR-SC',
    author_email='maximilian.betz@dlr.de',
    maintainer_email='martin.kuehn@dlr.de',
    url='https://github.com/SciCompMod/memilio',
    description='Part of MEmilio project, implementation of inference models for the existing models in MEmilio.',
    packages=find_packages(
        where=os.path.dirname(os.path.abspath(__file__))),
    install_requires=[
        # smaller pandas versions contain a bug that sometimes prevents reading
        # some excel files (e.g. population or twitter data)
        'pandas>=1.2.2',
        # smaller numpy versions cause a security issue, 1.25 breaks testing with pyfakefs
        'numpy>=1.22,<1.25',
        'bayesflow'
        'tensorflow',
        'matplotlib',
        'dataclasses',
        'dataclasses_json',],
    extras_require={'dev': [
        # first support of python 3.11
        'pyfakefs>=4.6',
        'coverage>=7.0.1',
    ], },
    long_description='', test_suite='memilio.inference_test',)
