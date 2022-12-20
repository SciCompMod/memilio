import os
import subprocess
import sys
from setuptools import setup, find_packages

try:
    from skbuild import setup
except ImportError:
    print('scikit-build is required to build from source.')
    print('Installation:  python -m pip install scikit-build')
    subprocess.check_call(
        [sys.executable, "-m", "pip", "install", "scikit-build"])
    from skbuild import setup

__version__ = '0.1.0'

setup(
    name='memilio-generation',
    version=__version__,
    author='DLR-SC',
    author_email='maximilian.betz@dlr.de',
    maintainer_email='martin.kuehn@dlr.de',
    url='https://github.com/DLR-SC/memilio',
    description='Part of MEmilio project, automatic generation of model specific python bindings.',
    packages=find_packages(
        where=os.path.dirname(os.path.abspath(__file__))),
    setup_requires=['cmake'],
    install_requires=['libclang', 'clang', 'dataclasses', 'dataclasses_json', ],
    extras_require={'dev': ['pyfakefs>=4.2.1', ], },
    long_description='', 
    test_suite='memilio.generation_test',
    package_data={'memilio': ['../_skbuild/linux-x86_64-3.8/cmake-build/compile_commands.json']},
)
