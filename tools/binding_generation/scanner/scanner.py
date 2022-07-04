import re
from os.path import exists

def get_infos(model, conf):

    # look if folder exists
    from subprocess import check_output
    project_path = check_output(["git", "rev-parse", "--show-toplevel"]).decode()[:-1] + "/cpp/models/"
    path = project_path + conf.folder

    # look if file exist
    try:
        f = open(path + "/model.h")
        f.close()
    except FileNotFoundError:
        print("File model.h not found. FileNotFoundError occured.")
    try:
        f = open(path + "/infection_state.h")
        f.close()
    except FileNotFoundError:
        print("File infection_state.h not found. FileNotFoundError occured.")
    try:
        f = open(path + "/parameters.h")
        f.close()
    except FileNotFoundError:
        print("File parameters.h not found. FileNotFoundError occured.")

    # search for infos
    namespaces = scan(path+"/model.h", r"^namespace", r"^namespace (?P<value>[a-z]+)", r"\{")
    model.namespace = (namespaces[0] +"::"+ namespaces[1] +"::")
    model.name = namespaces[1]

    population_groups = scan(path+"/model.h", r"class Model", r"Populations<(?P<value>[^>]+)>", r"\{")
    model.population_groups = population_groups[0].strip().split(",")

    model.init = scan(path+"/model.h", r"class (?P<model>[a-zA-Z]+) (?=: public CompartmentalModel)", r"Model\((?P<value>[^,]*)\)", r"void get_derivatives")

    model.compartments = scan(path+"/infection_state.h", r"enum class InfectionState", r"(?P<value>[a-zA-Z]+),", r"Count")




def scan(f_name, reg_str, data_reg_str, end_reg_str):

    reg = re.compile(reg_str)
    data_re = re.compile(data_reg_str)
    end_re = re.compile(end_reg_str)

    results = []

    with open(f_name) as file:
        need_search = False

        for line in file:
            if reg.search(line) is not None:
                need_search = True

            if need_search == True:
                res = data_re.search(line)
                if res is not None:
                    results.append(res.group('value'))

            if end_re.search(line) is not None:
                need_search = False

    return results