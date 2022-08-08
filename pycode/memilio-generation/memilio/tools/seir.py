from memilio.generation import Scanner, ScannerConfig, Generator
import os


here = os.path.dirname(os.path.abspath(__file__))
# Define ScannerConfig and initialize Scanner
with open(os.path.join(here + '/config.json'), 'r') as file:
    conf = ScannerConfig.schema().loads(file.read(), many=True)[0]
scanner = Scanner(conf)

# Extract results of Scanner into a intermed_repr
intermed_repr = scanner.extract_results()
#print(intermed_repr)

# Generate code
generator = Generator()
generator.create_substitutions(intermed_repr)
generator.generate_files(intermed_repr)

#scanner.output_ast_file()
