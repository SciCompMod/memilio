import sys
import os

BASE_PATH=os.path.dirname(__file__)

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
            'getrkidata=rki_data.getRKIData:main',
            'getpopuldata=rki_data.getPopulationData:main',
            'getjhdata = jh_data.getJHDataIntoPandasDataFrame:main',
            'getspaindata = data_spain.getSpainDataIntoPandasFrame:main'
        ],
    },
    package_dir={'': 'epidata'},
    packages=['rki_data', 'jh_data', 'data_spain'],
    long_description='',
    setup_requires=['cmake'],
    install_requires=requirements
)
