from generator.ModelManager import ModelManager
from generator.Model import Model
from generator.template.header import print_header
from generator.template.source import print_source

def parameter(name, c_type, default_value):
    """
    Create model parameter.

    Parameters
    ----------
    name - parameter name
    c_type - type of the parameter as string (e.g. "double")
    default_value - default parameter value

    Returns
    -------
    `sympy.Symbol` for this parameter
    """
    return ModelManager.get_active().add_parameter(name, c_type, default_value)

def compartments(*args):
    """
    Create model compartments and return corresponding symbols

    Parameters
    ----------
    *args : str or multiple str
        The names of model compartments, either space separated or as multiple strings

    Returns
    -------
    `sympy.Symbol` for each compartment
    """
    return ModelManager.get_active().set_compartments(*args)

def equation(lhs, op, rhs):
    """
    Create model equation. Equations are evaluated in order of creation.

    Parameters
    ----------
    lhs - `sympy.Symbol`
        symbol of a compartment or variable
    rhs - Expr
        any sympy expression
    op - string
        assignment operator, e.g. "=" or "+="
    """
    return ModelManager.get_active().add_equation(lhs, op, rhs)

def variables(*args):
    """
    Add one or more variables to the model and return corresponding symbols

    Parameters
    ----------
    *args : str or multiple str
        The names of variables, either space separated or as multiple strings

    Returns
    -------
    `sympy.Symbol` for each variable
    """
    return ModelManager.get_active().add_variable(*args)

def time():
    """
    creates the time variable to the model. only call once.

    Returns
    -------
    `sympy.Symbol` for time
    """
    # TODO: throw error when called twice
    return variables("t")

def model():
    """
    Create a new model. Usefull only if multiple models are created within the same scope. 

    Returns
    -------
    `generator.Model`
    """
    return ModelManager.add_model()

def delete(model):
    """
    Delete a model.

    Parameters
    ----------
    model : Model instance
    """
    err = str(model) + " is not of type Model"
    if type(model) is not Model: raise TypeError(err)
    ModelManager.remove_model(model)
    return

def name(string):
    """
    Set model name.

    Parameters
    ----------
    string - C compatible name
    """
    return ModelManager.get_active().set_name(string)

def namespace(string):
    """
    Set model namespace.

    Parameters
    ----------
    string - C compatible name
    """
    return ModelManager.get_active().set_namespace(string)
