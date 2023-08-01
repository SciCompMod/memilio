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
    install_requires=['libclang==14.0.6',
                      'dataclasses', 'dataclasses_json', 'importlib-resources>=1.1.0; python_version < \'3.9\''],
    extras_require={'dev': []},
    long_description='',
    test_suite='memilio.generation_test',
    package_data={'memilio.generation': [
        '../../_skbuild/*/cmake-build/compile_commands.json', '../tools/config.json']},
)
