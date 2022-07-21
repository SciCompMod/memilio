import generator as gen
from scanner.scanner_clang import Scanner
import scanner.scanner_config as config
import json

from scanner.scanner_config import ScannerConfig

# Define ScannerConfig and initialize Scanner
with open('config.json', 'r') as file:
    conf = ScannerConfig.schema().loads(file.read(), many=True)[0]
scanner = Scanner(conf)

# Extract results of Scanner into a model data
model = scanner.extract_results()
print(model)

# Generate code
generator = gen.Generator()
generator.create_substitutions(model)
generator.generate_files(model)

#scanner.output_ast_file()
