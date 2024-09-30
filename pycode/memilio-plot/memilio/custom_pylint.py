from setuptools import Command, Distribution
import os

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
            self.out_format)

    def run(self):
        os.makedirs("build_pylint", exist_ok=True)

        # Run pylint
        from pylint import lint
        with open(self.out_file, "w", encoding="utf-8") as report_file:
            options = ["--rcfile=../pylintrc", *self.lint_modules]

            lint.Run(options, reporter=self.reporter(
                report_file), do_exit=False)


def main():
    cmd = PylintCommand(dist=Distribution())
    cmd.initialize_options()
    cmd.finalize_options()
    cmd.run()


if __name__ == "__main__":
    main()
