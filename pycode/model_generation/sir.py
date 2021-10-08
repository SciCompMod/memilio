import generator as gen

project_path = ""

if project_path == "":
    from subprocess import check_output
    project_path = check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + "/cpp/"

s, i, r = gen.compartments("S I R")

beta = gen.parameter("ContactsPerDay", "double", 15)
gamma = gen.parameter("RecoveryRate", "double", 0.5)

gen.name("sir_example")
gen.namespace("Gen")

n = gen.variables("N")

gen.equation(n, "=", s + i + r)
gen.equation(s, "=", - beta/n * s * i)
gen.equation(i, "=", beta/n * s * i - gamma * i)
gen.equation(r, "=", gamma * i)

header_file = "epidemiology/secir/sir.h"
source_file = "examples/sir.cpp"

gen.print_header(project_path+header_file)
gen.print_source(header_file, 0, 1, 0.001, [9900, 100, 0], project_path+source_file)


