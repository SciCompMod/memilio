from setuptools import setup, find_packages

setup(
    name='memilio-stubs',
    version='0.1',
    packages=['memilio-stubs'],
    package_data={
        'memilio-stubs': ['simulation/*.pyi'],
    },
)
