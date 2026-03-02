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
        "name": "Parameters",
        "cursorkind": "CLASS_TEMPLATE",
        "methods": [
            "check_constraints",
            "apply_constraints",
            "get_start_commuter_detection",
            "get_end_commuter_detection",
            "get_commuter_nondetection"
        ]
    },

    {
        "type": "class",
        "name": "ModelBase",
        "cursorkind": "CLASS_TEMPLATE",
        "methods": []
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
    }


]
