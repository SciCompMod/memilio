from generator.utility import join
from generator.ModelManager import ModelManager

def print_source(headerfile, t0, tmax, dt, starting_population, file=None):
    m = ModelManager.get_active()
    m.finalize()
    s = source_string(m, headerfile, t0, tmax, dt, starting_population)
    if file is None:
        print(s)
    else:
        with open(file, 'w') as file:
            file.write(s)
    return

def source_string(model, headerfile, t0, tmax, dt, starting_population):
    assert(len(model.compartments) == len(starting_population))
    return (
        "{includes}\n"
        "{print_fct}\n"
        "int main() {{\n"
        "    epi::set_log_level(epi::LogLevel::debug);\n"
        "\n"
        "    // set start time, stop time and time step length\n"
        "    double t0   = {t0};\n"
        "    double tmax = {tmax};\n"
        "    double dt   = {dt};\n"
        "\n"
        "    epi::log_info(\"Simulating {model_name_caps}; t={{}} ... {{}} with dt = {{}}.\", t0, tmax, dt);\n"
        "\n"
        "    {namespace}::{model_name} model; // create {model_name} model instance\n"
        "\n"
        "    // set {model_name} model parameters"
        "    {parameter_settings}\n"
        "    // set starting population for each compartment"
        "    {population_settings}\n"
        "    printf(\"Starting population: {pop0_names}\\n\", {pop0_values});\n"
        "    printf(\"Total: %f\\n\", {pop0_sum});\n"
        "\n"
        "    // run {model_name} model simulation\n"
        "    epi::TimeSeries<double> {model_name} = epi::simulate(t0, tmax, dt, model);\n"
        "\n"
        "    print_to_terminal<double>({model_name}, std::vector<std::string>{{{compartments}}});\n"
        "}}\n"
    ).format(
        t0=t0, tmax=tmax, dt=dt,
        includes = include_string(model, headerfile),
        model_name      = model.name,
        model_name_caps = model.name.upper(),
        namespace = model.namespace.lower(),
        print_fct = print_function_string(),
        pop0_names =  join([c.name+" = %f" for c in model.compartments], ", "),
        pop0_values = join(["pop_"+c.name for c in model.compartments], ", "),
        pop0_sum = join(["pop_"+c.name for c in model.compartments], " + "),
        compartments = join(['"'+c.name+'"' for c in model.compartments], ", "),
        parameter_settings = parameter_settings_string(model),
        population_settings = populatioon_settings_string(model, starting_population)
    )

def parameter_settings_string(model, indent=4):
    tab = indent * " "
    s = "\r"
    for p in model.parameters:
        s += tab + (
                "model.parameters.set<{ns}::{pn}>({pd});\n"
            ).format(ns=model.namespace.lower(), pn=p.name, pd=p.default)
    return s

def populatioon_settings_string(model, starting_population, indent=4):
    tab = indent * " "
    s = "\r" + tab + "double "
    s += join(["pop_"+c.name+"="+str(p) for c, p in zip(model.compartments, starting_population)], ", ")
    s += ";\n"
    for c in model.compartments:
        s += tab + (
                "model.populations[{{epi::Index<{ns}::{pn}>({ns}::{pn}::{cn})}}] = pop_{cn};\n"
            ).format(ns=model.namespace.lower(), pn=model.compartment_name, cn=c.name)
    return s

def print_function_string():
    return (
        "#include <iostream>\n"
        "#include <vector>\n"
        "#include <string>\n"
        "\n"
        "template<class SC>\n"
        "void print_to_terminal(const epi::TimeSeries<SC>& results, const std::vector<std::string>& state_names) {{\n"
        "    printf(\"| %-16s |\", \"Time\");\n"
        "    for (size_t k = 0; k < state_names.size(); k++) {{\n"
        "        printf(\" %-16s |\", state_names[k].data()); // print underlying char*\n"
        "    }}\n"
        "    auto num_points = static_cast<size_t>(results.get_num_time_points());\n"
        "    for (size_t i = 0; i < num_points; i++) {{\n"
        "        printf(\"\\n| %16.6f |\", results.get_time(i));\n"
        "        auto res_i = results.get_value(i);\n"
        "        for (size_t j = 0; j < state_names.size(); j++) {{\n"
        "            printf(\" %16.6f |\", res_i[j]);\n"
        "        }}\n"
        "    }}\n"
        "    printf(\"\\n\");\n"
        "}}\n"
    ).format()

def include_string(model, headerfile):
    return (
        "#include \"{model_header}\"\n"
        "#include \"epidemiology/model/simulation.h\"\n"
    ).format(
        model_header = headerfile
    )