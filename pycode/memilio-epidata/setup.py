import os
import subprocess
import sys

from setuptools import Command, find_packages, setup

__version__ = '1.0.0'


class PylintCommand(Command):
    """
    Custom command to run pylint and get a report as html.
    """
    description = "Runs pylint and outputs the report as html."
    user_options = []

    def initialize_options(self):
        from pylint.reporters.json_reporter import JSONReporter
        from pylint.reporters.text import ParseableTextReporter, TextReporter
        from pylint_json2html import JsonExtendedReporter

        self.lint_modules = ["memilio/"]
        self.out_format = "extendedjson"

        self.REPORTERS = {
            "parseable": (ParseableTextReporter, "build_pylint/pylint_parseable.txt"),
            "text": (TextReporter, "build_pylint/pylint.txt"),
            "json": (JSONReporter, "build_pylint/pylint.json"),
            "extendedjson": (JsonExtendedReporter, "build_pylint/pylint_extended.json")
        }

    def finalize_options(self):
        self.reporter, self.out_file = self.REPORTERS.get(
            self.out_format)  # , self.REPORTERS.get("parseable"))

    def run(self):
        os.makedirs("build_pylint", exist_ok=True)

        # Run pylint
        from pylint import lint
        with open(self.out_file, "w", encoding="utf-8") as report_file:
            options = ["--rcfile=../pylintrc", *self.lint_modules]

            lint.Run(options, reporter=self.reporter(
                report_file), do_exit=False)


# Python-magic needs DLLs for libmagic. They have to be installed only on windows.
if sys.platform == 'win32':
    pymagic = 'python-magic-bin'
else:
    pymagic = 'python-magic'

setup(
    name='memilio-epidata',
    version=__version__,
    author='DLR-SC',
    author_email='daniel.abele@dlr.de',
    maintainer_email='martin.kuehn@dlr.de',
    url='https://github.com/SciCompMod/memilio',
    description='Part of MEmilio project, reads epidemiological data from different official and unofficial sources.',
    entry_points={
        'console_scripts': [
            'getcasedata=memilio.epidata.getCaseData:main',
            'getpopuldata=memilio.epidata.getPopulationData:main',
            'getjhdata = memilio.epidata.getJHData:main',
            'getdividata = memilio.epidata.getDIVIData:main',
            'getsimdata = memilio.epidata.getSimulationData:main',
            'cleandata = memilio.epidata.cleanData:main',
            'getcommutermobility = memilio.epidata.getCommuterMobility:main',
            'getvaccinationdata = memilio.epidata.getVaccinationData:main',
            'gethospitalizationdata = memilio.epidata.getHospitalizationData:main'
        ],
    },
    packages=find_packages(where=os.path.dirname(os.path.abspath(__file__))),
    long_description='',
    test_suite='memilio.epidata_test',
    install_requires=[
        # pandas 2.0 is minimum for CoW
        'pandas>=2.0.0',
        # FutureWarning of pandas that pyarrow will be required in a future release
        'pyarrow',
        'matplotlib',
        'tables',
        # smaller numpy versions cause a security issue, 1.25 breaks testing with pyfakefs
        'numpy>=1.22,<1.25',
        'openpyxl',
        'xlrd',
        'xlsxwriter',
        'requests',
        'pyxlsb',
        'wget',
        'twill==3.1',
        'python-calamine',
        pymagic
    ],
    extras_require={
        'dev': [
            # first support of python 3.11 4.6
            # 5.3.4 has conflicts with openpyxl
            # 5.3.3 broken
            'pyfakefs>=4.6,<5.3.3',
            # coverage 7.0.0 can't find .whl files and breaks CI
            'coverage>=7.0.1',
            # pylint 2.16 creates problem with wrapt package version
            'pylint>=2.13.0,<2.16',
            'pylint_json2html==0.4.0',
        ],
    },
    cmdclass={
        'pylint': PylintCommand
    },
)
