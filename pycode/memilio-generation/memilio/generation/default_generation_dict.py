default_dict = {
    "model": "Model",
    "agegroup": "AgeGroup",
    "emptystring": "",
    "simulation": "Simulation",
    "flowmodel": "FlowModel",
    "compartmentalmodel": "CompartmentaModel",
    "modelfile": "model.h",
    "infectionstatefile": "infection_state.h",
    "parameterspacefile": "parameter_space.h",
    "analyzeresultfile": "analyze_result.h",
    "namespace": "namespace",
    "mio": "mio",
    "drawsample": "draw_sample"
}

general_bindings_dict = [
    {
        "type": "class",
        "name": "Model",
        "cursorkind": "CLASS_TEMPLATE"
    },
    {
        "type": "function",
        "name": "draw_sample",
        "cursorkind": "FUNCTION_TEMPLATE",
    },
    {
        "type": "function",
        "name": "simulate_flows",
        "cursorkind": "FUNCTION_TEMPLATE"
    },
    {
        "type": "function",
        "name": "simulate",
        "cursorkind": "FUNCTION_TEMPLATE"
    },
    {
        "type": "function",
        "name": "check_constraints",
        "cursorkind": "CXX_METHOD"
    },
    {
        "type": "extern_function",
        "name": "interpolate_simulation_result",
        "namespace": "mio",
        "return_type": "mio::TimeSeries",
        "arg_types": ["const mio::TimeSeries<double>&", "const double"],
        "py_args": ["ts", "abs_tool"]
    },
    {
        "type": "extern_function",
        "name": "interpolate_simulation_result",
        "namespace": "mio",
        "return_type": "mio::TimeSeries",
        "arg_types": ["const mio::TimeSeries<double>&", "const std::vector<double>&"],
        "py_args": ["ts", "interpolation_times"]
    },

]
