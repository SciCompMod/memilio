#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
#
# Authors: Martin J. Kuehn
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy.stats import norm, gamma, lognorm, weibull_min, expon
plt.rc('legend', fontsize=14) 


# makes repeated lists unique, counts values and inserts zeros for values which were not there
def create_graph_from_repeated_list(repeated_list):
    unique, counts = np.unique(repeated_list, return_counts=True)
    
    if(len(unique)>0):   
        xvals = np.arange(max(unique)+1, dtype = int)
        yvals = np.zeros(max(unique)+1, dtype=int)

        for i in range(min(unique),max(unique)+1):
            yvals[unique] = counts

        if(yvals[0]==0):
            xvals = np.delete(xvals, 0) # account for distributions defined for x>0
            yvals = np.delete(yvals, 0)
            
    else:
        xvals = []
        yvals = []
        
    return (xvals, yvals)

    
# fits distribution to values given and checks via Pearsons Chi^2 test
# https://en.wikipedia.org/wiki/Pearson%27s_chi-squared_test
# n = the number of cells in the table
# The reduction in the degrees of freedom is calculated as p=s+1, 
# where s is the number of co-variates used in fitting the distribution (normal=2 (mu, sigma), lognormal=3...)
# returns argument 1&2: (x,y) to plot data,
#         argument 3&4: (x,y) to plot best fit, 
#         argument 5&6: distribution name and its parameters (shape (if existent), location and scale),
#         argument 7: SSE, chi2 score and pvalue in one row
def compute_distribution_fit(repeated_list, plot_fit=False):      
    args = []
    for i in range(len(repeated_list)):
        if repeated_list[i] < 0:
            args.append(i)
    repeated_list = np.delete(repeated_list, args)
    xvals, yvals = create_graph_from_repeated_list(repeated_list)
    yvals_scaled = 1/sum(yvals)*yvals
    
    xvals_refined = np.linspace(int(min(xvals)),max(repeated_list)+1,200)     
    
    params_list = []
    
    #https://docs.scipy.org/doc/scipy/reference/stats.html#continuous-distributions
    dist_names = ['gamma', 'weibull_min', 'lognorm', 'norm', 'expon']
    error_list = np.zeros([len(dist_names),3])

    i=0
    for dist_name in dist_names:
        # get name of distribution
        dist_name_scipy = getattr(scipy.stats, dist_name)

        # fit dist to data
        params = dist_name_scipy.fit(repeated_list)
        
        # separate parts of parameters
        shape = params[:-2] # shape parameter for e.g. gamma, lognorm and weibull distributions
                            # parameter not available for normal and exponential
        loc = params[-2] # location parameter
        scale = params[-1] # scale parameter
        
        params_list.append(params)        
        
        y_fitted = dist_name_scipy.pdf(xvals, *shape, loc=loc, scale=scale)
        sse = np.power(yvals_scaled - y_fitted, 2.0) # compute squared errors
        
        # compute chi2 and p-value
        chi2 = scipy.stats.chisquare(yvals_scaled, f_exp=y_fitted) 
        sse = np.sum(sse) # compute sum 
        
        error_list[i,0] = sse
        error_list[i,1] = chi2[0] # chi2 score
        error_list[i,2] = chi2[1] # pvalue
        
        i += 1
        
        
    best_fit_dist_ind = np.argsort(error_list[:,1], axis=0)
    
    for i in range(0,1):
        # separate parts of parameters
            # shape parameter for e.g. gamma, lognorm and weibull distributions
            # parameter not available for normal and exponential
        shape = params_list[best_fit_dist_ind[i]][:-2] 
        loc = params_list[best_fit_dist_ind[i]][-2] # location parameter
        scale = params_list[best_fit_dist_ind[i]][-1] # scale parameter

        dist_name_scipy = getattr(scipy.stats, dist_names[best_fit_dist_ind[i]])

        y_fitted_plot = dist_name_scipy.pdf(xvals_refined, *shape, loc=loc, scale=scale)
        
        if plot_fit:
            label_plot = str(dist_names[best_fit_dist_ind[i]])+\
            '('+str(np.around(*shape,2))+','+str(np.around(loc,2))+','+str(np.around(scale,2))+')'

            plt.plot(xvals_refined, y_fitted_plot, label=label_plot)
            plt.xlim(min(repeated_list),max(repeated_list))
        
    
    if plot_fit:
        plt.plot(xvals, yvals_scaled, label='original data')
        plt.legend(loc='upper right')
        plt.show()  
        print('Best fitted distribution: ',dist_names[best_fit_dist_ind[0]])
        print('\tSSE: ',error_list[best_fit_dist_ind[0],0], \
              "\n\tChi2: ",error_list[best_fit_dist_ind[0],1], \
              "\n\tpvalue: ",error_list[best_fit_dist_ind[0],2])
    
    return (xvals, yvals, xvals_refined, y_fitted_plot, \
            dist_names[best_fit_dist_ind[0]], \
            params_list[best_fit_dist_ind[i]], \
            error_list[best_fit_dist_ind[0],:])