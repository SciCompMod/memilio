import generator as gen

project_path = ""

if project_path == "":
    from subprocess import check_output
    project_path = check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + "/pycode/memilio-simulation/memilio/simulation/"


gen.name("oseir")
gen.namespace("mio::oseir::")

header_file = "oseir.py"
source_file = "oseir.cpp"

gen.print_python_file(project_path+header_file)
gen.print_source(project_path+source_file)