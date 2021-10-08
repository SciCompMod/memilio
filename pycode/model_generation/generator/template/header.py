
from sympy import symbols as Symbols
from sympy.printing import ccode
from generator.utility import join, assign
from generator.ModelManager import ModelManager

def print_header(file=None):
    m = ModelManager.get_active()
    m.finalize()
    s = header_string(m)
    if file is None:
        print(s)
    else:
        with open(file, 'w') as file:
            file.write(s)
    return

def header_string(model):
    return (
        "#ifndef {model_name_caps}_H_\n"
        "#define {model_name_caps}_H_\n"
        "\n"
        "{includes}\n"
        "namespace {namespace} {{\n"
        "\n"
        "{population}\n"
        "{parameterset}\n"
        "class {model_name} : public epi::CompartmentalModel<{population_name}, {parameterset_name}> {{\n"
        "public:\n"
        "    {model_name}() : epi::CompartmentalModel<{population_name}, {parameterset_name}>("
        "{population_name}(epi::Index<{compartment_name}>({compartment_name}::Count)), "
        "{parameterset_name}()) {{}};\n"
        "\n"
        "    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop,\n"
        "                         Eigen::Ref<const Eigen::VectorXd> y, double t,\n"
        "                         Eigen::Ref<Eigen::VectorXd> dydt) const override\n"
        "    {{\n"
        "        const {parameterset_name}& par = this->parameters;\n"
        "        {computation}"
        "    }};\n"
        "}};\n"
        "\n"
        "}} // end namespace {namespace}\n"
        "\n"
        "#endif"
    ).format(
        model_name_caps   = model.namespace.upper() + "_" + model.name.upper(),
        namespace         = model.namespace.lower(),
        model_name        = model.name,
        includes          = include_string(),
        population        = population_string(model.compartments, model.compartment_name, model.population_name),
        population_name   = model.population_name,
        compartment_name  = model.compartment_name,
        parameterset      = parameterset_string(model.parameters, model.parameterset_name),
        parameterset_name = model.parameterset_name,
        computation       = computation_string(model, model.compartment_name, 8)
    )

def parameterset_string(parameters, parameterset_name):
    if len(parameters) == 0:
        return ""
    s = ""
    p_names = join([p.name for p in parameters], ", ")
    for p in parameters:
        s += (
                "struct {struct_name} {{\n"
                "    using Type = {dtype};\n"
                "    static Type get_default() {{\n"
                "        return {default_value};\n"
                "    }}\n"
                "}};\n"
            ).format(dtype=p.c_type, struct_name=p.name, default_value=p.default)
    s += ("\nusing {n} = epi::ParameterSet<{p}>;\n").format(p=p_names, n=parameterset_name)
    return s

def population_string(compartment_list, compartment_name, population_name):
    s = compartment_list[0].name
    for c in compartment_list[1:]:
        s += ",\n    " + c.name
    return (
        "enum {n} {{\n"
        "    {comps},\n"
        "    Count\n"
        "}};\n"
        "\n"
        "using {p} = epi::Populations<{n}>;\n"
    ).format(comps=s, n=compartment_name, p=population_name)

def include_string():
    return (
        '#include "epidemiology/math/smoother.h"\n'
        '#include "epidemiology/model/compartmentalmodel.h"\n'
        '#include "epidemiology/model/populations.h"\n'
        '#include "epidemiology/model/simulation.h"\n'
        '#include "epidemiology/utils/parameter_set.h"\n'
        '#include "epidemiology/utils/eigen_util.h"\n'
        '#include "epidemiology/utils/index.h"\n'
    )

def computation_string(model, compartment_name, indent):
    tab = indent * " "
    s = "\r"
    c_lhs = [Symbols("dydt["+compartment_name+r"\:\:"+c.name+"]") for c in model.compartments]
    c_rhs = [Symbols("pop["+compartment_name+r"\:\:"+c.name+"]") for c in model.compartments]
    param = [Symbols("par.get<"+p.name+">()") for p in model.parameters]
    var_decl = [Symbols(r"auto\ " + v.name) for v in model.variables]

    for e in model.equations:
        for c, c_dydt, c_pop in zip(model.compartments, c_lhs, c_rhs):
            e = assign(e.lhs.subs(c.symbol, c_dydt), e.rhs.subs(c.symbol, c_pop), e.op)
        for p, par in zip(model.parameters, param):
            e = e.subs(p.symbol, par)
        for i in range(len(model.variables)):
            v = model.variables[i]
            if (v.symbol == e.lhs) and (var_decl[i] != None):
                e = assign(e.lhs.subs(v.symbol, var_decl[i]), e.rhs, e.op)
                var_decl[i] = None
        # set human=False to avoid warnings in the printed code (e.h. "//Not supported in c")
        s += tab + ccode(e, human=False)[2] + "\n"
    return s