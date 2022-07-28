from memilio.generation import Scanner, ScannerConfig, Generator
import json

# Define ScannerConfig and initialize Scanner
with open('./memilio/tools/config.json', 'r') as file:
    conf = ScannerConfig.schema().loads(file.read(), many=True)[0]
#scanner = Scanner(conf)

# Extract results of Scanner into a model data
#model = scanner.extract_results()
#print(model)

# Generate code
#generator = Generator()
#generator.create_substitutions(model)
#generator.generate_files(model)

#scanner.output_ast_file()
