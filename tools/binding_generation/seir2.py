import generator as gen
import scanner.scanner as scanner
import scanner.scanner_config as config

project_path = ""

if project_path == "":
    from subprocess import check_output
    project_path = check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + "/pycode/memilio-simulation/memilio/simulation/"

conf = config.ScannerConfig()
conf.folder = "ode_seir"
conf.model_name = "test_seir"

model = gen.Model()
scanner.get_infos(model, conf)

model.finalize(conf)
print(model)

header_file = model.pymio_name + ".py"
source_file = model.pymio_name + ".cpp"

gen.print_python_file(model, project_path+header_file)
gen.print_source(model, project_path+source_file)