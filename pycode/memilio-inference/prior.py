import numpy as np
from scipy import stats
from enum import Enum, auto
from functools import partial

from bayesflow.simulation import Prior

alpha_f = (0.7**2)*((1-0.7)/(0.17**2) - (1-0.7))
beta_f = alpha_f*(1/0.7 - 1)

# Global variable to hold the ParamNames enum
ParamNames = None


class ModelStrategy:
    """Abstract base class for different model strategies."""

    @staticmethod
    def generate_base():
        raise NotImplementedError

    @staticmethod
    def generate_intervention():
        raise NotImplementedError

    @staticmethod
    def generate_observation(sim_diff):
        raise NotImplementedError


class SIRStrategy(ModelStrategy):
    param_names_base = [
        r"$\lambda_0$", r"$\mu$", r"$I_0$"]
    param_names_intervention = [
        # r"$\Delta t_1$", r"$\Delta t_2$", r"$\Delta t_3$", r"$\Delta t_4$",
        r"$t_1$", r"$t_2$", r"$t_3$", r"$t_4$",
        r"$\lambda_1$", r"$\lambda_2$", r"$\lambda_3$", r"$\lambda_4$"]
    param_names_observation = [
        r"$f_i$", r"$\phi_i$", r"$D_i$", r"$\psi$"]

    @staticmethod
    def generate_base():  # Possible option to draw without redraw
        """Generates a random draw from the joint prior."""

        def check_if_parameter_constraints_satisfied(params):
            # Extract parameters
            lambd0, mu, I0 = params

            # Impose constraints
            try:
                assert 1 >= lambd0 >= 0
                # Check for mu
            except AssertionError as e:
                return False
            return True

        while True:  # somehow check that we dont get into softlock
            lambd0 = np.random.lognormal(mean=np.log(1.2), sigma=0.5)
            # mu = np.random.lognormal(mean=np.log(1/8), sigma=0.2) # distribution of paper, need 1/dist
            mu = np.random.lognormal(mean=np.log(8), sigma=0.2)
            I0 = np.random.gamma(shape=2, scale=30)

            # return draws if they are inside domains and follow certain constraints
            if check_if_parameter_constraints_satisfied([lambd0, mu, I0]):
                return [lambd0, mu, I0]

    @staticmethod
    def generate_intervention():
        """Generates a random draw from the joint prior."""

        def check_if_parameter_constraints_satisfied(params):
            # Extract parameters
            t1, t2, t3, t4, lambd1, lambd2, lambd3, lambd4 = params
            t1, t2, t3, t4 = int(round(t1)), int(round(t2)), int(
                round(t3)), int(round(t4))  # is that changing the params?

            # Impose constraints
            try:
                assert t1 > 0 and t2 > 0 and t3 > 0 and t4 > 0
                # do i really want to make this constraint, or just test for last t, because its a special case
                assert t1 < t2 < t3 < t4
                # assert delta_t1 > 0 and delta_t2 > 0 and delta_t3 > 0 and delta_t4 > 0
                # assert t2 - t1 >= delta_t1 and t3 - t2 >= delta_t2 and t4-t3 >= delta_t3 and T-t4 >= delta_t4 # i dont use delta_t for now
                assert 1 >= lambd1 >= 0 and 1 >= lambd2 >= 0 and 1 >= lambd3 >= 0 and 1 >= lambd4 >= 0
            except AssertionError as e:
                return False
            return True

        while True:  # somehow check that we dont get into softlock
            t1 = np.random.normal(loc=8, scale=3)
            t2 = np.random.normal(loc=15, scale=1)
            t3 = np.random.normal(loc=22, scale=1)
            t4 = np.random.normal(loc=66, scale=1)
            # delta_t1 = np.random.lognormal(mean=np.log(3), sigma=0.3)
            # delta_t2 = np.random.lognormal(mean=np.log(3), sigma=0.3)
            # delta_t3 = np.random.lognormal(mean=np.log(3), sigma=0.3)
            # delta_t4 = np.random.lognormal(mean=np.log(3), sigma=0.3)
            lambd1 = np.random.lognormal(mean=np.log(0.6), sigma=0.5)
            lambd2 = np.random.lognormal(mean=np.log(0.3), sigma=0.5)
            lambd3 = np.random.lognormal(mean=np.log(0.1), sigma=0.5)
            lambd4 = np.random.lognormal(mean=np.log(0.1), sigma=0.5)

            # return draws if they are inside domains and follow certain constraints
            if check_if_parameter_constraints_satisfied([t1, t2, t3, t4,  # delta_t1, delta_t2, delta_t3, delta_t4,
                                                        lambd1, lambd2, lambd3, lambd4]):
                return [t1, t2, t3, t4,  # delta_t1, delta_t2, delta_t3, delta_t4,
                        lambd1, lambd2, lambd3, lambd4]

    @staticmethod
    def generate_observation(sim_diff):
        """Generates a random draw from the joint prior."""

        def check_if_parameter_constraints_satisfied(params):
            # Extract parameters
            f_i, phi_i, D_i, scale_I = params

            # Impose constraints
            try:
                # need simdiff as context
                assert (sim_diff - 1) > int(round(D_i))
                # Check for f_i, phi_i, scale_I
                assert scale_I > 0
                assert 1 >= f_i >= 0
            except AssertionError as e:
                return False
            return True

        while True:  # somehow check that we dont get into softlock
            f_i = np.random.beta(a=alpha_f, b=beta_f)
            phi_i = stats.vonmises(kappa=0.01).rvs()
            D_i = np.random.lognormal(mean=np.log(8), sigma=0.2)
            scale_I = np.random.gamma(shape=1, scale=5)

            # return draws if they are inside domains and follow certain constraints
            if check_if_parameter_constraints_satisfied([f_i, phi_i, D_i, scale_I]):
                return [f_i, phi_i, D_i, scale_I]


def set_param_enum(param_names):
    """Function to create the param Enum dynamically"""
    global ParamNames  # Use the global keyword to modify the global variable
    param_dict = {
        name.replace(r"$", "").replace("\\", "").replace("{", "").replace("}", "").replace("_", ""): auto()
        for name in param_names
    }
    ParamNames = Enum('ParamNames', param_dict)


# def generate_builder(sim_diff=16):
#     """
#     First draft of a builder for model prior. For now just single option without different builders,
#     extend into different possibilities or combine with data_generation model.
#     """
#     param_names = [
#         r"$\lambda_0$", r"$\mu$", r"$I_0$",
#         # r"$\Delta t_1$", r"$\Delta t_2$", r"$\Delta t_3$", r"$\Delta t_4$",
#         r"$t_1$", r"$t_2$", r"$t_3$", r"$t_4$",
#         r"$\lambda_1$", r"$\lambda_2$", r"$\lambda_3$", r"$\lambda_4$",
#         r"$f_i$", r"$\phi_i$", r"$D_i$", r"$\psi$"
#     ]

#     def model_prior():
#         params = []
#         params.extend(generate_base())
#         params.extend(generate_intervention())
#         params.extend(generate_observation(sim_diff))
#         return params
#     return model_prior, param_names


# def fixed_model_params_wrapper(fixed_params=dict()):
#     """
#     First draft of giving an easy interface to fix certain parameters.
#     In future maybe combine this into a big model prior builder.
#     """
#     def decorator_function(original_function):
#         def wrapper_function(*args, **kwargs):
#             params = original_function(*args, **kwargs)
#             for param_name, value in fixed_params.items():
#                 # If the key is a string, convert it to the corresponding enum value
#                 if isinstance(param_name, str):
#                     param_enum = ParamNames[param_name]
#                 elif isinstance(param_name, ParamNames):
#                     param_enum = param_name
#                 else:
#                     raise ValueError(f"Invalid parameter name: {param_name}")

#                 index = param_enum.value - 1  # Get the index of the parameter
#                 params[index] = value
#             return params
#         return wrapper_function
#     return decorator_function


class ModelPriorBuilder:
    def __init__(self, strategy, sim_diff=16):
        self.strategy = strategy
        self.sim_diff = sim_diff
        self.components = []
        self.param_names = []
        self.fixed_parameters = {}

    def add_base(self):
        self.components.append(self.strategy.generate_base)
        self.param_names.extend(self.strategy.param_names_base)
        set_param_enum(self.param_names)
        return self

    def add_intervention(self):
        self.components.append(self.strategy.generate_intervention)
        self.param_names.extend(self.strategy.param_names_intervention)
        set_param_enum(self.param_names)
        return self

    def add_observation(self):
        self.components.append(
            partial(self.strategy.generate_observation, self.sim_diff))
        self.param_names.extend(self.strategy.param_names_observation)
        set_param_enum(self.param_names)
        return self

    def set_fixed_parameters(self, fixed_params=None):
        self.fixed_parameters = fixed_params
        return self

    def build(self):
        def model_prior():
            params = []
            for component in self.components:
                params.extend(component())
            for param_name, value in self.fixed_parameters.items():
                # If the key is a string, convert it to the corresponding enum value
                if isinstance(param_name, str):
                    param_enum = ParamNames[param_name]
                elif isinstance(param_name, ParamNames):
                    param_enum = param_name
                else:
                    raise ValueError(f"Invalid parameter name: {param_name}")

                index = param_enum.value - 1
                params[index] = value
            return params
        return Prior(prior_fun=model_prior, param_names=self.param_names)


if __name__ == "__main__":
    # Building a SIR model
    sir_builder = ModelPriorBuilder(SIRStrategy())
    prior = (sir_builder.add_base()
                        .add_intervention()
                        .add_observation()
                        .set_fixed_parameters({ParamNames.I0: 0, ParamNames.mu: 0.5})
                        .build())

    params_sir = prior(batch_size=1)
    print(params_sir)
