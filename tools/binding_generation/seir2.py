import generator as gen
from scanner.scanner_clang import Scanner
import scanner.scanner_config as config

project_path = ""

if project_path == "":
    from subprocess import check_output
    project_path = check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + "/pycode/memilio-simulation/memilio/simulation/"

# Define ScannerConfig and initialize Scanner
conf = config.ScannerConfig()
conf.set_folder("ode_seir")
conf.set_python_module_name("test_seir")
scanner = Scanner(conf)

# Extract results of Scanner into a model data
model = gen.Model()
scanner.extract_results(model)
model.finalize(conf)
#print(model)

# Generate code
generator = gen.Generator()
generator.create_substitutions(model)
generator.generate_files(model)

scanner.output_ast(2)
