import os
from setuptools import setup, find_packages

__version__ = '0.1.0'

setup(
    name='memilio-generation',
    version=__version__,
    author='DLR-SC',
    author_email='maximilian.betz@dlr.de',
    maintainer_email='daniel.abele@dlr.de',
    url='https://github.com/DLR-SC/memilio',
    description='Part of MEmilio project, python bindings to the C++ libraries that contain the models and simulations.',
    packages=find_packages(where=os.path.dirname(os.path.abspath(__file__))),
    install_requires=[
        'sympy',
        'matplotlib',
        'libclang',
        'clang',
        'dataclasses_json',
    ],
    extras_require={
        'dev': [
            'numpy >= 1.21',
            'pyfakefs>=4.2.1',
        ],
    },
    long_description='',
    test_suite='memilio.generation_test',
)
