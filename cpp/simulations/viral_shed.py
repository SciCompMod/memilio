import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
fontsize = 24

def viral_load(t, increase, decrease, peak):
    peak_time = peak / increase
    if t < peak_time:
        return increase * t
    else:
        return peak - decrease * (t - peak_time)
    
class viral_load_with_params:
    def __init__(self, increase, decrease, peak):
        self.increase = increase
        self.decrease = decrease
        self.peak = peak
        
    def __call__(self, t):
        return viral_load(t, self.increase, self.decrease, self.peak)

def viral_shedding(t, personal_viral_shed, alpha, beta, viral_load_with_params):
    return personal_viral_shed * 1/(1+np.exp(-(alpha + beta* viral_load_with_params(t))))

def return_my_sigma(mean, stddev):
    my = np.log(mean * mean / np.sqrt(mean * mean + stddev * stddev));
    sigma = np.sqrt(np.log(1 + stddev * stddev / mean / mean));
    return my, sigma

def draw_timepoints_infectious_and_not():

    mean_incubation = 4.5
    stddev_incubation = 1.5
    my_sigma_incubation = return_my_sigma(mean_incubation, stddev_incubation)
    exposed_length = sp.stats.lognorm(s=my_sigma_incubation[0], scale=my_sigma_incubation[1]).rvs(1)



    
    return timepoints_infectious, timepoints_not_infectious

viral_load_increase = 2.0
viral_load_decrease = 0.5
viral_load_peak = 8.1
personal_viral_shed = sp.stats.gamma(a=1.6, scale=1/22).rvs(20)
alpha = -7
beta = 1

viral_load_with_params = viral_load_with_params(viral_load_increase, viral_load_decrease, viral_load_peak)

t = np.linspace(0.0000001, 20, 100)
viral_shedding_values_for_each_person = [[viral_shedding(t, personal_viral_shed, alpha, beta, viral_load_with_params) for t in t] for personal_viral_shed in personal_viral_shed]

plt.figure(figsize=(12, 8))
for viral_shedding_values in viral_shedding_values_for_each_person:
    plt.plot(t, viral_shedding_values)
plt.title('Viral shedding', fontsize=fontsize)
plt.xlabel('Time (days)', fontsize=fontsize-4)
plt.ylabel('Viral shed', fontsize=fontsize-4)
plt.setp(plt.gca().get_xticklabels(), fontsize=fontsize-8)  
plt.setp(plt.gca().get_yticklabels(), fontsize=fontsize-8)

# we need to save the figure with high resolution and jpeg format
plt.savefig('/Users/saschakorf/Downloads/viral_shedding.jpg', dpi=300)


# we do the same for viral load
viral_load_values = [np.exp(viral_load_with_params(t)) for t in t]

plt.figure(figsize=(12, 8))
plt.plot(t, viral_load_values)
plt.title('Viral load', fontsize=fontsize)
plt.xlabel('Time (days)', fontsize=fontsize-4)
plt.ylabel('Viral load', fontsize=fontsize-4)
plt.setp(plt.gca().get_xticklabels(), fontsize=fontsize-8)
plt.setp(plt.gca().get_yticklabels(), fontsize=fontsize-8)

# we need to save the figure with high resolution and jpeg format
plt.savefig('viral_load.jpg', dpi=300)

