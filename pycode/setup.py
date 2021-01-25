import sys
import os
import subprocess
import distutils.cmd

try:
    from skbuild import setup
except ImportError:
    print('scikit-build is required to build from source.')
    print('Installation:  python -m pip install scikit-build')
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scikit-build"])
    from skbuild import setup


__version__ = '0.1.0'


INSTALL_REQUIRES = ['pandas<1.2.0',
                    'matplotlib',
                    'tables',
                    'numpy<=1.19.4',
                    'openpyxl',
                    'xlrd']

EXTRAS_REQUIRE = {"pylint": ["pylint", "pylint_json2html"]}

class PylintCommand(distutils.cmd.Command):
    """
    Custom command to run pylint and get a report as html.
    """
    description = "Runs pylint and outputs the report as html."
    user_options = []

    def initialize_options(self):
        from pylint.reporters.text import TextReporter, ParseableTextReporter
        from pylint.reporters.json_reporter import JSONReporter
        from pylint_json2html import JsonExtendedReporter

        self.lint_modules = ["epidemiology/"]

        self.out_format = "extendedjson"

        self.REPORTERS = {
            "parseable": (ParseableTextReporter, "build_pylint/pylint_parseable.txt"),
            "text": (TextReporter, "build_pylint/pylint.txt"),
            "json": (JSONReporter, "build_pylint/pylint.json"),
            "extendedjson": (JsonExtendedReporter, "build_pylint/pylint_extended.json")
        }

    def finalize_options(self):
        self.reporter, self.out_file = self.REPORTERS.get(self.out_format)#, self.REPORTERS.get("parseable"))

    def run(self):
        os.makedirs("build_pylint", exist_ok=True)

        # Run pylint
        from pylint import lint
        with open(self.out_file, "w", encoding="utf-8") as report_file:
            options = ["--rcfile=pylintrc", "-j 2", *self.lint_modules]

            lint.Run(options, reporter=self.reporter(report_file), do_exit=False)

setup(
    name='epidemiology',
    version=__version__,
    author='DLR-SC',
    author_email='martin.siggel@dlr.de',
    maintainer_email='kathrin.rack@dlr.de',
    url='https://gitlab.dlr.de/hpc-against-corona/epidemiology',
    description='The python package for the HPC corona project',
    entry_points={
        'console_scripts': [
            'getrkidata=epidemiology.epidata.getRKIData:main',
            'getpopuldata=epidemiology.epidata.getPopulationData:main',
            'getjhdata = epidemiology.epidata.getJHData:main',
            'getspaindata = epidemiology.epidata.getSpainData:main',
            'getdividata = epidemiology.epidata.getDIVIData:main',
            'getalldata = epidemiology.epidata.getAllData:main',
            'cleandata = epidemiology.epidata.cleanData:main',
            'getrkiestimation = epidemiology.epidata.getRKIDatawithEstimations:main'
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
    test_suite='test',
    install_requires=INSTALL_REQUIRES,
    extras_require = EXTRAS_REQUIRE,
    cmdclass={
            'pylint': PylintCommand,
        },
)
