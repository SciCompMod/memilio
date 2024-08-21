import numpy as np
from scipy.stats import lognorm


def get_lognormal_parameters(mean, std):
    """
    Compute shape and scale parameters to use in lognormal distribution for given mean and standard deviation. 
    The lognormal distribution we consider in state_age_function.h is based on the implementation in scipy and the parameters 
    shape and scale are defined accordingly. 
    """
    variance = std**2

    mean_tmp = np.log(mean**2/np.sqrt(mean**2+variance))
    variance_tmp = np.log(variance/mean**2 + 1)

    shape = np.sqrt(variance_tmp)
    scale = np.exp(mean_tmp)

    # Test if mean and std are as expected for computed shape and scale parameters.
    mean_lognorm, variance_lognorm = lognorm.stats(
        shape, loc=0, scale=scale, moments='mv')

    mean_test = np.exp(scale**2/2)

    # print("Test mean: ", mean_test)

    # print(mean_lognorm, variance_lognorm)

    if np.abs(mean_lognorm-mean) > 1e-8:
        print('Distribution does not have expected mean value.')

    if np.abs(np.sqrt(variance_lognorm)-std) > 1e-8:
        print('Distribution does not have expected standard deviation.')

    return round(shape, 8), round(scale, 8)


def get_weighted_mean(prob_1, stay_time_1, stay_time_2):

    weighted_mean = prob_1*stay_time_1 + (1-prob_1)*stay_time_2

    return weighted_mean


if __name__ == '__main__':
    shape, scale = get_lognormal_parameters(2.183, 1.052)
    print(f"{shape:.12f}", f"{scale:.12f}")

    weighted_mean = get_weighted_mean(0.793099, 1.1, 8.0)
    print(f"{weighted_mean:.6f}")

    weighted_mean = get_weighted_mean(0.078643, 6.6, 8.0)
    print(f"{weighted_mean:.6f}")

    weighted_mean = get_weighted_mean(0.173176, 1.5, 18.1)
    print(f"{weighted_mean:.6f}")

    weighted_mean = get_weighted_mean(0.387803, 10.7, 18.1)
    print(f"{weighted_mean:.6f}")
