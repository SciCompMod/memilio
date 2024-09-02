import operator
from typing import Self, Any
from collections.abc import Callable
from inspect import isclass

import numpy as np
from bayesflow.simulation import Prior


class PriorScaler():

    def fit(self: Self, prior: Prior) -> Self:
        prior_means, prior_stds = prior.estimate_means_and_stds()
        for idx, prior_std in enumerate(prior_stds[0]):
            if prior_std == 0:
                prior_stds[0, idx] += 10 ** -14
        self.prior_means = prior_means
        self.prior_stds = prior_stds
        return self

    def transform(self: Self, data: Any) -> Any:
        # z-standardize with previously computed means and stds
        # Make this use different scalers
        self._check_is_fitted()
        return (data - self.prior_means) / self.prior_stds

    def inverse_transform(self: Self, data: Any) -> Any:
        self._check_is_fitted()
        return self.prior_means + data * self.prior_stds

    def _check_is_fitted(self: Self):
        if isclass(self):
            raise TypeError(
                f"{self} is a class, not an instance.")
        msg = (
            "This %(name)s instance is not fitted yet. Call 'fit' with "
            "appropriate arguments before using this estimator."
        )
        for attribute in ["prior_means", "prior_stds"]:
            if not hasattr(self, attribute):
                raise AttributeError(msg % {"name": type(self).__name__})


class UnboundParameter():
    def __init__(self: Self, distribution: Callable[[], float], name: str) -> None:
        self.distribution = distribution
        self.name = name

    def get_draw(self: Self) -> Callable[[], float]:
        def draw() -> float:
            return self.distribution()
        return draw


class FixedParameter(UnboundParameter):
    def __init__(self: Self, value, name: str) -> None:
        super().__init__(value, name)

    def get_draw(self: Self) -> Callable[[], float]:
        def draw() -> float:
            return self.distribution
        return draw


class BoundParameter(UnboundParameter):
    def __init__(self: Self, distribution: Callable[[], float], name: str, lower_bound: float = None,
                 upper_bound: float = None, lower_bound_included: bool = True, upper_bound_included: bool = True) -> None:
        super().__init__(distribution, name)

        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.lower_bound_operator = operator.le if lower_bound_included else operator.lt
        self.upper_bound_operator = operator.ge if upper_bound_included else operator.gt

    def check(self: Self, value) -> bool:
        if self.lower_bound is not None and not self.lower_bound_operator(self.lower_bound, value):
            return False
        if self.upper_bound is not None and not self.upper_bound_operator(self.upper_bound, value):
            return False
        return True

    def get_draw(self: Self) -> Callable[[], float]:
        def draw():
            while True:  # somehow check that we dont get into softlock
                value = self.distribution()
                if self.check(value):
                    return value
        return draw


class LambdaParameter(BoundParameter):
    def __init__(self: Self, distribution: Callable[[], float], name: str) -> None:
        # check 1 >= lambd >= 0
        super().__init__(distribution, name, lower_bound=0, upper_bound=1,
                         lower_bound_included=True, upper_bound_included=True)


class InterventionChangePointParameter(BoundParameter):
    def __init__(self: Self, distribution: Callable[[], float], name: str) -> None:
        # ckeck change_point >= 0
        super().__init__(distribution, name, lower_bound=0, lower_bound_included=True)


class DelayParameter(BoundParameter):
    def __init__(self: Self, distribution: Callable[[], float], name: str, sim_diff: int) -> None:
        # check ((sim_diff - 1) > int(round(delay))), thus -1.5
        super().__init__(distribution, name,
                         upper_bound=sim_diff - 1.5, upper_bound_included=False)


class WeeklyModulationAmplitudeParameter(BoundParameter):
    def __init__(self: Self, distribution: Callable[[], float], name: str) -> None:
        # check 1 >= f_i >= 0
        super().__init__(distribution, name, lower_bound=0, upper_bound=1,
                         lower_bound_included=True, upper_bound_included=True)


class ScaleMultiplicativeReportingNoiseParameter(BoundParameter):
    def __init__(self: Self, distribution: Callable[[], float], name: str) -> None:
        # check scale_I > 0
        super().__init__(distribution, name, lower_bound=0, lower_bound_included=False)


class ModelStrategy:
    """Abstract base class for different model strategies."""

    @staticmethod
    def add_base(prior_array: list[UnboundParameter]) -> None:
        raise NotImplementedError

    @staticmethod
    def add_intervention(prior_array: list[UnboundParameter]) -> None:
        raise NotImplementedError

    @staticmethod
    def add_observation(prior_array: list[UnboundParameter], sim_diff: int) -> None:
        raise NotImplementedError


class ModelPriorBuilder:
    def __init__(self: Self, strategy: ModelStrategy, sim_diff: int = 16):
        self.strategy = strategy
        self.sim_diff = sim_diff
        self.components = []
        self.fixed_parameters = []

    def add_base(self: Self) -> Self:
        self.strategy.add_base(self.components)
        return self

    def add_intervention(self: Self) -> Self:
        self.strategy.add_intervention(self.components)
        return self

    def add_observation(self: Self) -> Self:
        self.strategy.add_observation(self.components, self.sim_diff)
        return self

    def add_dummy(self: Self) -> Self:
        self.components.append(UnboundParameter(
            distribution=np.random.uniform, name=r'Dummy'))
        return self

    def set_fixed_parameters(self: Self, fixed_params: dict[str, float] = None) -> Self:
        self.fixed_parameters = []  # maybe delete fixed params and just save names
        for param_name, value in fixed_params.items():
            self.fixed_parameters.append(
                FixedParameter(value=value, name=param_name))
        return self

    def build(self: Self) -> Prior:
        parameter_priors = []
        param_names = []

        # get components
        for component in self.components:
            #   skip fixed params
            if any(component.name == fixed_param.name for fixed_param in self.fixed_parameters):
                continue
            parameter_priors.append(component.get_draw())
            param_names.append(component.name)

        # # redefine components with scalar value
        # for fixed_param in self.fixed_parameters:
        #     param_name = fixed_param.name
        #     idx = param_names.index(param_name)
        #     parameter_priors[idx] = fixed_param.get_draw()

        def model_prior() -> list[float]:
            parameter_draws = []
            for prior in parameter_priors:
                parameter_draws.append(prior())
            return parameter_draws
        return Prior(prior_fun=model_prior, param_names=param_names)


# class SIRStrategy(ModelStrategy):
#     param_names_base = [
#         r"$\lambda_0$", r"$\mu$", r"$I_0$"]
#     param_names_intervention = [
#         # r"$\Delta t_1$", r"$\Delta t_2$", r"$\Delta t_3$", r"$\Delta t_4$",
#         r"$t_1$", r"$t_2$", r"$t_3$", r"$t_4$",
#         r"$\lambda_1$", r"$\lambda_2$", r"$\lambda_3$", r"$\lambda_4$"]
#     param_names_observation = [
#         r"$f_i$", r"$\phi_i$", r"$D_i$", r"$\psi$"]

#     @staticmethod
#     def generate_base():  # Possible option to draw without redraw
#         """Generates a random draw from the joint prior."""

#         def check_if_parameter_constraints_satisfied(params):
#             # Extract parameters
#             lambd0, mu, I0 = params

#             # Impose constraints
#             try:
#                 assert 1 >= lambd0 >= 0
#                 # Check for mu
#             except AssertionError as e:
#                 return False
#             return True

#         while True:  # somehow check that we dont get into softlock
#             lambd0 = np.random.lognormal(mean=np.log(1.2), sigma=0.5)
#             # mu = np.random.lognormal(mean=np.log(1/8), sigma=0.2) # distribution of paper, need 1/dist
#             mu = np.random.lognormal(mean=np.log(8), sigma=0.2)
#             I0 = np.random.gamma(shape=2, scale=30)

#             # return draws if they are inside domains and follow certain constraints
#             if check_if_parameter_constraints_satisfied([lambd0, mu, I0]):
#                 return [lambd0, mu, I0]

#     @staticmethod
#     def generate_intervention():
#         """Generates a random draw from the joint prior."""

#         def check_if_parameter_constraints_satisfied(params):
#             # Extract parameters
#             t1, t2, t3, t4, lambd1, lambd2, lambd3, lambd4 = params
#             t1, t2, t3, t4 = int(round(t1)), int(round(t2)), int(
#                 round(t3)), int(round(t4))  # is that changing the params?

#             # Impose constraints
#             try:
#                 assert t1 > 0 and t2 > 0 and t3 > 0 and t4 > 0
#                 # do i really want to make this constraint, or just test for last t, because its a special case
#                 assert t1 < t2 < t3 < t4
#                 # assert delta_t1 > 0 and delta_t2 > 0 and delta_t3 > 0 and delta_t4 > 0
#                 # assert t2 - t1 >= delta_t1 and t3 - t2 >= delta_t2 and t4-t3 >= delta_t3 and T-t4 >= delta_t4 # i dont use delta_t for now
#                 assert 1 >= lambd1 >= 0 and 1 >= lambd2 >= 0 and 1 >= lambd3 >= 0 and 1 >= lambd4 >= 0
#             except AssertionError as e:
#                 return False
#             return True

#         while True:  # somehow check that we dont get into softlock
#             t1 = np.random.normal(loc=8, scale=3)
#             t2 = np.random.normal(loc=15, scale=1)
#             t3 = np.random.normal(loc=22, scale=1)
#             t4 = np.random.normal(loc=66, scale=1)
#             # delta_t1 = np.random.lognormal(mean=np.log(3), sigma=0.3)
#             # delta_t2 = np.random.lognormal(mean=np.log(3), sigma=0.3)
#             # delta_t3 = np.random.lognormal(mean=np.log(3), sigma=0.3)
#             # delta_t4 = np.random.lognormal(mean=np.log(3), sigma=0.3)
#             lambd1 = np.random.lognormal(mean=np.log(0.6), sigma=0.5)
#             lambd2 = np.random.lognormal(mean=np.log(0.3), sigma=0.5)
#             lambd3 = np.random.lognormal(mean=np.log(0.1), sigma=0.5)
#             lambd4 = np.random.lognormal(mean=np.log(0.1), sigma=0.5)

#             # return draws if they are inside domains and follow certain constraints
#             if check_if_parameter_constraints_satisfied([t1, t2, t3, t4,  # delta_t1, delta_t2, delta_t3, delta_t4,
#                                                         lambd1, lambd2, lambd3, lambd4]):
#                 return [t1, t2, t3, t4,  # delta_t1, delta_t2, delta_t3, delta_t4,
#                         lambd1, lambd2, lambd3, lambd4]

#     @staticmethod
#     def generate_observation(sim_diff):
#         """Generates a random draw from the joint prior."""

#         def check_if_parameter_constraints_satisfied(params):
#             # Extract parameters
#             f_i, phi_i, D_i, scale_I = params

#             # Impose constraints
#             try:
#                 # need simdiff as context
#                 assert (sim_diff - 1) > int(round(D_i))
#                 # Check for f_i, phi_i, scale_I
#                 assert scale_I > 0
#                 assert 1 >= f_i >= 0
#             except AssertionError as e:
#                 return False
#             return True

#         while True:  # somehow check that we dont get into softlock
#             f_i = np.random.beta(a=alpha_f, b=beta_f)
#             phi_i = stats.vonmises(kappa=0.01).rvs()
#             D_i = np.random.lognormal(mean=np.log(8), sigma=0.2)
#             scale_I = np.random.gamma(shape=1, scale=5)

#             # return draws if they are inside domains and follow certain constraints
#             if check_if_parameter_constraints_satisfied([f_i, phi_i, D_i, scale_I]):
#                 return [f_i, phi_i, D_i, scale_I]


# class ModelPriorBuilder:
#     def __init__(self, strategy, sim_diff=16):
#         self.strategy = strategy
#         self.sim_diff = sim_diff
#         self.components = []
#         self.param_names = []
#         self.fixed_parameters = {}

#     def add_base(self):
#         self.components.append(self.strategy.generate_base)
#         self.param_names.extend(self.strategy.param_names_base)
#         # set_param_enum(self.param_names)
#         return self

#     def add_intervention(self):
#         self.components.append(self.strategy.generate_intervention)
#         self.param_names.extend(self.strategy.param_names_intervention)
#         # set_param_enum(self.param_names)
#         return self

#     def add_observation(self):
#         self.components.append(
#             partial(self.strategy.generate_observation, self.sim_diff))
#         self.param_names.extend(self.strategy.param_names_observation)
#         # set_param_enum(self.param_names)
#         return self

#     def set_fixed_parameters(self, fixed_params=None):
#         self.fixed_parameters = fixed_params
#         return self

#     def build(self):
#         def model_prior():
#             params = []
#             for component in self.components:
#                 params.extend(component())
#             for param_name, value in self.fixed_parameters.items():
#                 # If the key is a string, convert it to the corresponding enum value
#                 if isinstance(param_name, str):
#                     index = self.param_names.index(param_name)
#                 else:
#                     raise ValueError(f"Invalid parameter name: {param_name}")
#                 params[index] = value
#             return params
#         return Prior(prior_fun=model_prior, param_names=self.param_names)
