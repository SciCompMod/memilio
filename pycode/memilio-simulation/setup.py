import os

from setuptools import find_packages

from skbuild import setup

setup(
    packages=find_packages(where=os.path.dirname(os.path.abspath(__file__))),
    setup_requires=['cmake'],
    # need shared libs so there is one shared log level
    cmake_args=['-DMEMILIO_BUILD_SHARED_LIBS:BOOL=ON'],
    install_requires=[
        # smaller pandas versions contain a bug that sometimes prevents reading
        # some excel files (e.g. population or twitter data)
        'pandas>=2.0.0',
    ],
    extras_require={
        'dev': [
            # smaller numpy versions cause a security issue, 1.25 breaks testing with pyfakefs
            'numpy>=1.22,<1.25',
        ],
    },
    long_description='', test_suite='memilio.simulation_test',
)
