import os
import subprocess
import sys

from setuptools import find_packages, setup

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
    name='memilio-simulation',
    version=__version__,
    author='DLR-SC',
    author_email='daniel.abele@dlr.de',
    maintainer_email='daniel.abele@dlr.de',
    url='https://github.com/DLR-SC/memilio',
    description='Part of MEmilio project, python bindings to the C++ libraries that contain the models and simulations.',
    packages=find_packages(where=os.path.dirname(os.path.abspath(__file__))),
    setup_requires=['cmake'],
    install_requires=[],
    extras_require={
        'dev': ['numpy >= 1.22'],
    },
    long_description='',
    test_suite='memilio.simulation_test',
)
