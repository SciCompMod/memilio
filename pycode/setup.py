import sys
import os
import subprocess


try:
    from skbuild import setup
except ImportError:
    print('scikit-build is required to build from source.')
    print('Installation:  python -m pip install scikit-build')
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scikit-build"])
    from skbuild import setup


__version__ = '0.1.0'


INSTALL_REQUIRES = ['pandas', 'matplotlib', 'tables']


setup(
    name='epidemiology',
    version=__version__,
    author='DLR-SC',
    author_email='martin.siggel@dlr.de',
    url='https://gitlab.dlr.de/hpc-against-corona/epidemiology',
    description='The python package for the HPC corona project',
    entry_points={
        'console_scripts': [
            'getrkidata=epidemiology.epidata.getRKIData:main',
            'getpopuldata=epidemiology.epidata.getPopulationData:main',
            'getjhdata = epidemiology.epidata.getJHData:main',
            'getspaindata = epidemiology.epidata.getSpainData:main'
        ],
    },
    package_dir={
       'epidemiology': 'epidemiology',
       'epidemiology.seir': os.path.join('epidemiology', 'seir'),
       'epidemiology.secir': os.path.join('epidemiology','secir'),
       'epidemiology.epidata': os.path.join('epidemiology','epidata')},
    packages=['epidemiology', 
              'epidemiology.seir',
              'epidemiology.secir',
              'epidemiology.epidata'],
    long_description='',
    setup_requires=['cmake'],
    install_requires=INSTALL_REQUIRES
)
