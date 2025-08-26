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

__version__ = '1.0.0'

setup(
    name='memilio-simulation', version=__version__, author='DLR-SC',
    author_email='daniel.abele@dlr.de', maintainer_email='Martin.Kuehn@DLR.de',
    url='https://github.com/SciCompMod/memilio',
    description='Part of MEmilio project, python bindings to the C++ libraries that contain the models and simulations.',
    packages=find_packages(where=os.path.dirname(os.path.abspath(__file__))),
    # setup_requires=['cmake'],
    # need shared libs so there is one shared log level
    cmake_args=['-DMEMILIO_BUILD_SHARED_LIBS:BOOL=ON'],
    install_requires=[
    ],
    extras_require={
        'dev': [
            # smaller numpy versions cause a security issue, 1.25 breaks testing with pyfakefs
            'numpy>=1.22,<1.25',
            # smaller pandas versions contain a bug that sometimes prevents reading
            # some excel files (e.g. population or mobility data)
            'pandas>=2.0.0',
        ],
    },
    long_description='', test_suite='memilio.simulation_test',
)
