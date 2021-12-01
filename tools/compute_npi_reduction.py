import numpy as np
import os

# path and read only works if repo is called 'memilio'
data_path = os.getcwd().split('memilio')[0] + 'memilio/data/'
print("My current directory is : " + data_path)

# names for contact location data files
home = 'home'
school = 'school_pf_eig'
work = 'work'
other = 'other'

# list of names
contact_locs = [home, school, work, other]

# if minimum contact pattern > 0 is used; if false, minimum = 0
use_minimum = False

# provide population of age groups as used in contact matrices
population = np.array([3961376,7429883,19117865,28919134,18057318,5681135])

# read contact matrices
cm_bl = []
cm_bl_min = []
for loc in contact_locs:
    cm_bl.append(np.loadtxt(data_path + '/contacts/baseline_' + loc + '.txt'))
    if use_minimum:
        cm_bl_min.append(np.loadtxt(data_path + '/contacts/minimum_' + loc + '.txt'))
    else:
        cm_bl_min.append(np.zeros((6,6)))
    
# compute total number of contacts
contacts_total_base = []
contacts_total_min = []
for i in range(len(contact_locs)):
    contacts_total_base.append(sum(np.matmul(cm_bl[i], population)))
    contacts_total_min.append(sum(np.matmul(cm_bl_min[i], population)))

# print out basic info
print('Number of baseline in home: '
      + str(np.round(100*contacts_total_base[0] / sum(contacts_total_base), 2)) + '%, school: '
      + str(np.round(100*contacts_total_base[1] / sum(contacts_total_base), 2)) + '%, work: '
      + str(np.round(100*contacts_total_base[2] / sum(contacts_total_base), 2)) + '%, other: '
      + str(np.round(100*contacts_total_base[3] / sum(contacts_total_base), 2)) + '%')

# define intervals of contact reduction for different levels
reduc_factors_min = np.array([[0.5, 0.3, 0.6, 0.6], [0.0, 0.25, 0.25, 0.25]])
reduc_factors_max = np.array([[0.7, 0.5, 0.8, 0.8], [0.0, 0.35, 0.35, 0.35]])

# make three copies and compute for min, mean and max
contacts_total_reduced = [contacts_total_base.copy(
), contacts_total_base.copy(), contacts_total_base.copy()]
if reduc_factors_min.shape == reduc_factors_max.shape:
    for i in range(len(contacts_total_base)):
        reduc = [1, 1, 1]
        for j in range(reduc_factors_min.shape[0]):
            # compute combined (all levels) reduction factors
            reduc[0] *= (1 - reduc_factors_min[j][i]) # min
            reduc[1] *= (1 - 0.5*(reduc_factors_min[j][i] + reduc_factors_max[j][i])) # mean
            reduc[2] *= (1 - reduc_factors_max[j][i]) # max
        for k in range(len(reduc)):
            contacts_total_reduced[k][i] = contacts_total_reduced[k][i] - (1 - reduc[k]) * (contacts_total_reduced[k][i] - contacts_total_min[i])

    print('Summed contact reduction (including protection effects): ' + \
          str(np.round(100*(1 - sum(contacts_total_reduced[1]) / sum(contacts_total_base)), 2))
          + '% (' + str(np.round(100*(1 - sum(contacts_total_reduced[0]) / sum(contacts_total_base)), 2)) 
          + ' - ' + str(np.round(100*(1 - sum(contacts_total_reduced[2]) / sum(contacts_total_base)), 2)) + '%)')

    print('Average number of contact reduction in home: '
        + str(np.round(100*(1 - contacts_total_reduced[1][0] / contacts_total_base[0]), 2)) + '%, school: '
        + str(np.round(100*(1 - contacts_total_reduced[1][1] / contacts_total_base[1]), 2)) + '%, work: '
        + str(np.round(100*(1 - contacts_total_reduced[1][2] / contacts_total_base[2]), 2)) + '%, other: '
        + str(np.round(100*(1 - contacts_total_reduced[1][3] / contacts_total_base[3]), 2)) + '%')

    print('Number of unprotected contacts after reduction in home: '
        + str(np.round(100*contacts_total_reduced[1][0] / sum(contacts_total_reduced[1]), 2)) + '%, school: '
        + str(np.round(100*contacts_total_reduced[1][1] / sum(contacts_total_reduced[1]), 2)) + '%, work: '
        + str(np.round(100*contacts_total_reduced[1][2] / sum(contacts_total_reduced[1]), 2)) + '%, other: '
        + str(np.round(100*contacts_total_reduced[1][3] / sum(contacts_total_reduced[1]), 2)) + '%')

else:
    print('Error in reduction factor dimension. Please provide equally '
        'sized arrays for minimum and maximum.')

