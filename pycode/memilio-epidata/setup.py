import sys
import os
import subprocess
from setuptools import setup, find_packages, Command

__version__ = '0.1.0'


class PylintCommand(Command):
    """
    Custom command to run pylint and get a report as html.
    """
    description = "Runs pylint and outputs the report as html."
    user_options = []

    def initialize_options(self):
        from pylint.reporters.text import TextReporter, ParseableTextReporter
        from pylint.reporters.json_reporter import JSONReporter
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


setup(
    name='memilio-epidata',
    version=__version__,
    author='DLR-SC',
    author_email='daniel.abele@dlr.de',
    maintainer_email='martin.kuehn@dlr.de',
    url='https://github.com/DLR-SC/memilio',
    description='Part of MEmilio project, reads epidemiological data from different official and unofficial sources.',
    entry_points={
        'console_scripts': [
            'getcasedata=memilio.epidata.getCaseData:main',
            'getpopuldata=memilio.epidata.getPopulationData:main',
            'getjhdata = memilio.epidata.getJHData:main',
            'getdividata = memilio.epidata.getDIVIData:main',
            'getsimdata = memilio.epidata.getSimulationData:main',
            'cleandata = memilio.epidata.cleanData:main',
            'getcasesestimation = memilio.epidata.getCaseDatawithEstimations:main',
            'getcommutermobility = memilio.epidata.getCommuterMobility:main'
        ],
    },
    packages=find_packages(where=os.path.dirname(os.path.abspath(__file__))),
    long_description='',
    test_suite='memilio.epidata_test',
    install_requires=[
        'pandas<1.2.0',
        'matplotlib<3.4',
        'tables',
        'numpy>=1.21',
        'openpyxl',
        'xlrd',
        'requests',
	    'pyxlsb',
        'wget'
    ],
    extras_require={
        'dev': [
            'pyfakefs==4.1.0',
            'freezegun',
            'coverage',
            'pylint<=2.11.1',
            'pylint_json2html<=0.3.0',
        ],
    },
    cmdclass={
        'pylint': PylintCommand
    },
)
