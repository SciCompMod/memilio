import sys

try:
    from skbuild import setup
except ImportError:
    print('scikit-build is required to build from source.', file=sys.stderr)
    print('Please run:', file=sys.stderr)
    print('', file=sys.stderr)
    print('  python -m pip install scikit-build')
    sys.exit(1)

__version__ = '0.1.0'

requirements = ['pandas', 'matplotlib', 'tables']

setup(
    name='epidemiology',
    version=__version__,
    author='DLR-SC',
    author_email='martin.siggel@dlr.de',
    url='https://gitlab.dlr.de/hpc-against-corona/epidemiology',
    description='The python package for the HPC corona project',
    entry_points={
        'console_scripts': [
            'getrkidata=epidemiology.epidata.getGeoJsonIntoPandasDataFrame_for_RKI:cli',
        ],
    },
    package_dir={'epidemiology.epidata': 'epidata'},
    packages=['epidemiology.epidata'],
    long_description='',
    setup_requires=['cmake'],
    install_requires=requirements
)
