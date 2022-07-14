import generator as gen
from scanner.scanner_clang import Scanner
import scanner.scanner_config as config
import json

project_path = ""

if project_path == "":
    from subprocess import check_output
    project_path = check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + "/pycode/memilio-simulation/memilio/simulation/"

# Define ScannerConfig and initialize Scanner
file = open("config.json")
conf = json.load(file)
scanner = Scanner(conf)

# Extract results of Scanner into a model data
#model = scanner.extract_results()
#model.finalize(conf)
#print(model)

# Generate code
#generator = gen.Generator()
#generator.create_substitutions(model)
#generator.generate_files(model)

scanner.output_ast_file()
