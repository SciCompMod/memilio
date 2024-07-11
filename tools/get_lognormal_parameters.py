import numpy as np
from scipy.stats import lognorm


def get_lognormal_parameters(mean, std):
    """
    Compute shape and scale parameters to use in lognormal distribution for given mean and standard deviation. 
    The lognormal distribution we consider in state_age_function.h is based on the implementation in scipy and the parameters 
    shape and scale are defined accordingly. 
    """
    variance = std**2

    mean_tmp = np.log(mean) - 0.5*np.log(1 + variance/mean**2)
    variance_tmp = np.log(variance/mean**2 + 1)

    shape = np.sqrt(variance_tmp)
    scale = np.exp(mean_tmp)

    # Test if mean and std are as expected for computed shape and scale parameters.
    mean_lognorm, variance_lognorm = lognorm.stats(
        shape, loc=0, scale=scale, moments='mv')

    if np.abs(mean_lognorm-mean) > 1e-8:
        print('Distribution does not have expected mean value.')

    if np.abs(np.sqrt(variance_lognorm)-std) > 1e-8:
        print('Distribution does not have expected standard deviation.')

    return shape, scale


def get_weighted_mean(prob_1, stay_time_1, stay_time_2):

    weighted_mean = prob_1*stay_time_1 + (1-prob_1)*stay_time_2

    return weighted_mean


if __name__ == '__main__':
    shape, scale = get_lognormal_parameters(10.7, 4.8)
    print(f"{shape:.8f}", f"{scale:.8f}")

    weighted_mean = get_weighted_mean(0.217177, 10.7, 18.1)
    print(f"{weighted_mean:.8f}")
