
import random 
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D



def plot_dist(all_dampings, name):
    value_count = pd.Series(np.asarray(all_dampings).flatten()).value_counts().sort_index()
    days_list= [*range(0,len(value_count)+1, 1)] 
    plt.clf()
    plt.stairs(value_count, days_list, fill = True)
    #plt.hist(value_count)
    plt.show()
    plt.title('Distribution Random Sampling')
    plt.xlabel('Days')
    plt.ylabel('Count')
    plt.savefig("random"+name+".png")


###########################################
# neuer Ansatz mit Martin Schattendampings 
# best method --> not perfektly uniform but when we perform a limited number of runs ( e.g 10k) the results are acceptable 
# this method rather punishes the first and last few days
# in our case this is better than having more dampings towards the end of the prediction horizon (as it is slighly the case for the circle method)

def generate_dampings_withshadowdamp():

    number_of_dampings = 5
    days = 90
    min_distance = 7 
    min_damping_day = 0
    number_of_runs= 10000

    all_dampings = []
    count_runs = 0 
    count_shadow = 0
    while len(all_dampings)<number_of_runs:

        days_list = list(range((min_damping_day),days))
        dampings = []
        if count_shadow <2:
            for i in range(number_of_dampings):
                damp = random.choice(days_list)
                days_before = list(range(damp-(min_distance), damp))
                days_after = list(range(damp, damp+(min_distance+1)))
                dampings.append(damp)
                days_list = [ele for ele in days_list if ele not in (days_before+days_after)]
        else: 
            # chose a forbidden damping 
            damp = random.choice(list(range((0-min_distance),0))+ list(range(days+1, (days+min_distance+1))))
                
            days_before = list(range(damp-(min_distance), damp))
            days_after = list(range(damp, damp+(min_distance+1)))
            days_list = [ele for ele in days_list if ele not in (days_before+days_after)]
            dampings.append(damp)
            for i in range(number_of_dampings):
                
                
                damp = random.choice(days_list)
                days_before = list(range(damp-(min_distance), damp))
                days_after = list(range(damp, damp+(min_distance+1)))
                dampings.append(damp)
                days_list = [ele for ele in days_list if ele not in (days_before+days_after)]
                count_shadow = 0
        
            
        forbidden_damping_values = list(range((0-min_distance),0))+ list(range(days+1, (days+min_distance+1)))
        dampings = [ele for ele in dampings if ele not in forbidden_damping_values]
        count_runs+=1
        count_shadow +=1
        # select first or last five dampings
        if len(dampings) >= number_of_dampings:
            #dampings = random.sample(dampings, 5)
            all_dampings.append(sorted(dampings))
        #     if count_runs % 2 == 0:
        
        return all_dampings



#### for comparison: without shadow damps

def normal_drawing():
    number_of_dampings = 2
    days = 100
    min_distance = 7 
    min_damping_day = 5
    number_of_runs= 500000

    all_dampings = []
    while len(all_dampings)<number_of_runs:
        count_runs = 0 
        days_list = list(range((min_damping_day-min_distance),days+(min_distance+1)))
        dampings = []
        for i in range(number_of_dampings):
            damp = random.choice(days_list)
            days_before = list(range(damp-(min_distance+1), damp))
            days_after = list(range(damp, damp+(min_distance+1)))
            dampings.append(damp)
            days_list = [ele for ele in days_list if ele not in (days_before+days_after)]
            
        forbidden_damping_values = list(range((min_damping_day-min_distance),min_damping_day))+ list(range(days+1, (days+min_distance+1)))
        dampings = [ele for ele in dampings if ele not in forbidden_damping_values]
        count_runs+=1
        
        # select first or last five dampings
        if len(dampings) >= number_of_dampings:
            all_dampings.append(sorted(dampings))
        #     if count_runs % 2 == 0:
                #all_dampings.append(sorted(dampings[:number_of_dampings]))
        #     else: 
        #         all_dampings.append(sorted(dampings[-number_of_dampings:]))





##### circle method

replace_dict = {-7:84, -6:85, -5:86, -4:87, -3:88, -2:89, -1:90, 91:0, 92:1, 93:2, 94:3, 95:4, 96:5, 97:6}

number_of_dampings = 5
days = 90
min_distance = 7 
min_damping_day = 1
number_of_runs= 1000000

all_dampings = []
count_runs = 0 

while len(all_dampings)<number_of_runs:

    days_list = list(range(min_damping_day-1,days))
    dampings = []
 
    for i in range(number_of_dampings):
            damp = random.choice(days_list)
            days_before = list(range(damp-(min_distance), damp))
            days_before_replaced =  [(replace_dict[d] if d in replace_dict else d) for d in days_before]
            days_after = list(range(damp, damp+(min_distance+1)))
            days_after_replaced =  [(replace_dict[d] if d in replace_dict else d) for d in days_after]

            dampings.append(damp)
            days_list = [ele for ele in days_list if ele not in (days_before_replaced+days_after_replaced)]
    
            

    #forbidden_damping_values = min_damping_day-1
    #forbidden_damping_values = list(range((0-min_distance),0))+ list(range(days+1, (days+min_distance+1)))
    #dampings = [ele for ele in dampings if ele not in forbidden_damping_values]
    count_runs+=1
   
   
    if len(dampings) >= 5:
        #dampings = random.sample(dampings, 5)
        all_dampings.append(sorted(dampings))






##########################  plots  ##########################################################
# plot average number of days between dampings 
def plot_differences(all_dampings, name):
    differences = []
    for days_array in all_dampings: 
        differences.append(np.diff(days_array).tolist())

    plt.clf()
    #plt.stairs(value_count, days_list, fill = True)
    plt.hist(np.asarray(differences).flatten(), bins = len(np.unique(np.asarray(differences).flatten())))
    plt.show()
    plt.title('Distribution of Distances to next Damping')
    plt.xlabel('Days')
    plt.ylabel('Count')
    plt.savefig("random_differences" + name+ ".png")



#### boxplot per damping 
def plot_boxplot(all_dampings, name):
    plt.clf()
    # Set the figure size
    plt.rcParams["figure.figsize"] = [7.50, 3.50]
    plt.rcParams["figure.autolayout"] = True

    data = pd.DataFrame(data = all_dampings, columns = [1,2,3,4,5])

    # Plot the dataframe
    ax = data.plot(kind='box')

    plt.title('Boxplots for each Damping')
    plt.xlabel('Damping Number')
    plt.ylabel('Damping Day')
    plt.savefig("random_boxplots"+name+ ".png")

# get stats 

#print(data_info)

def plot_maxminmean(all_dampings, name):
    data = pd.DataFrame(data = all_dampings, columns = [1,2,3,4,5])
    d = {'Mean': data.mean().values.astype(int), 'Min': data.min().values, 'Max':data.max().values}
    data_info = pd.DataFrame(data = d)
    df = data_info
    plt.figure().clf() 

    linestyles = ['--', '-', ':']

    markers = []
    for m in Line2D.markers:
        try:
            if len(m) == 1 and m != " ":
                markers.append(m)
        except TypeError:
            pass

    for i, ls, m in zip(df.columns, linestyles, markers): 
        plt.plot(df[i], label = i,  linestyle=ls, marker = m )


    plt.title('Mean, Max and Min per Damping')
    plt.legend()
    plt.xlabel('Damping Number')
    plt.xticks(data.columns.values)
    plt.ylabel('Days')
    plt.savefig("maxminmean" + name + ".png")



# plot average number of days between dampings 
def plot_differences(all_dampings, name):
    differences = []
    for days_array in all_dampings: 
        differences.append(np.diff(days_array).tolist())

    plt.clf()
    #plt.stairs(value_count, days_list, fill = True)
    plt.hist(np.asarray(differences).flatten(), bins = len(np.unique(np.asarray(differences).flatten())))
    plt.show()
    plt.title('Distribution of Distances to next Damping')
    plt.xlabel('Days')
    plt.ylabel('Count')
    plt.savefig("random_differences" + name+ ".png")



#### boxplot per damping 
def plot_boxplot(all_dampings, name):
    plt.clf()
    # Set the figure size
    plt.rcParams["figure.figsize"] = [7.50, 3.50]
    plt.rcParams["figure.autolayout"] = True

    data = pd.DataFrame(data = all_dampings, columns = [1,2,3,4,5])

    # Plot the dataframe
    ax = data.plot(kind='box')

    plt.title('Boxplots for each Damping')
    plt.xlabel('Damping Number')
    plt.ylabel('Damping Day')
    plt.savefig("random_boxplots"+name+ ".png")

# get stats 

#print(data_info)

def plot_maxminmean(all_dampings, name):
    data = pd.DataFrame(data = all_dampings, columns = [1,2,3,4,5])
    d = {'Mean': data.mean().values.astype(int), 'Min': data.min().values, 'Max':data.max().values}
    data_info = pd.DataFrame(data = d)
    df = data_info
    plt.figure().clf() 

    linestyles = ['--', '-', ':']

    markers = []
    for m in Line2D.markers:
        try:
            if len(m) == 1 and m != " ":
                markers.append(m)
        except TypeError:
            pass

    for i, ls, m in zip(df.columns, linestyles, markers): 
        plt.plot(df[i], label = i,  linestyle=ls, marker = m )


    plt.title('Mean, Max and Min per Damping')
    plt.legend()
    plt.xlabel('Damping Number')
    plt.xticks(data.columns.values)
    plt.ylabel('Days')
    plt.savefig("maxminmean" + name + ".png")

name = ' final_'
plot_dist(all_dampings, days, name )




######################################### alte Ansätze ####################################
all_dampings = []
number_of_dampings = 5
puffer_damp = 2
runs = 9000000
min_damping_day = 0 
days = 90
min_distance = 7 
for i in range (runs):
    damping_days=(random.sample(range(min_damping_day-min_distance,days+min_distance), number_of_dampings+puffer_damp))
    all_dampings.append(sorted(damping_days))



#calculate difference
differences = []
indices = []
for days_array, index in zip(all_dampings, np.arange(0, len(all_dampings)) ): 
    differences.append(np.diff(days_array).tolist())
    indices.append(index)

number_of_differences = len(differences[0])

res = [[i, d] for [i,d] in zip(indices, differences)  if len([j for j in d if j >= min_distance]) == number_of_differences]



df_diff = pd.DataFrame(data = res, columns = ['index', 'diff'])

cleaned_all_dampings = []
for index in sorted(df_diff['index'], reverse=True):
    cleaned_all_dampings.append(all_dampings[index])



forbidden_damping_values = list(range((0-min_distance),0))+ list(range(days+1, (days+min_distance+1)))


cleaned_all_dampings_withoutforbidden = [] 
for i in cleaned_all_dampings :
        cleaned_all_dampings_withoutforbidden.append([ele for ele in i  if ele not in forbidden_damping_values])



###########################################################################

days = 90
number_of_runs = runs/10
normal_damp = []
while len(normal_damp)<number_of_runs:
    count_runs = 0 
    days_list = list(range((0),days))
    dampings = []
    for i in range(number_of_dampings):
        damp = random.choice(days_list)
        days_before = list(range(damp-(min_distance+1), damp))
        days_after = list(range(damp, damp+(min_distance+1)))
        dampings.append(damp)
        days_list = [ele for ele in days_list if ele not in (days_before+days_after)]
        
    count_runs+=1
    
    # select first or last five dampings
    if len(dampings) >= 5:
        normal_damp.append(sorted(dampings))


final = []
for i in cleaned_all_dampings_withoutforbidden:
       if len(i) >= 5:
            dampings =  random.sample(i, k = 5 )
            final.append(sorted(dampings))


final = final + normal_damp


################ completely random, only preventing duplicates #################
all_dampings = []
for i in range (10):
    damping_days=(random.sample(range(0,days), number_of_dampings))
    all_dampings.append(sorted(damping_days))




########### generate randomly and keep what fulfills condition ########

# generate more than enough data, since we know, that only about 20% of the random generator is suitable to our conditions 
all_dampings_2 = []
runs = 8000000
for i in range (runs):
    damping_days=(random.sample(range(min_damping_day,days), number_of_dampings))
    all_dampings_2.append(sorted(damping_days))

# delete all samples where at least one damping is less than 7 days apart from the other one 

#calculate difference
differences = []
indices = []
for days_array, index in zip(all_dampings_2, np.arange(0, len(all_dampings_2)) ): 
    differences.append(np.diff(days_array).tolist())
    indices.append(index)

number_of_differences = len(differences[0])

res = [[i, d] for [i,d] in zip(indices, differences)  if len([j for j in d if j >= min_distance]) == number_of_differences]

#res = [[i, d] for [i,d] in zip(indices, differences)  if len([j for j in d if j >= min_distance]) == number_of_differences]
df_diff = pd.DataFrame(data = res, columns = ['index', 'diff'])

cleaned_all_dampings = []
for index in sorted(df_diff['index'], reverse=True):
    cleaned_all_dampings.append(all_dampings_2[index])


###################################################################################################################

# add for (min_distance)/(days+min_distance) of all cases a damping in the days days to days+min_distance and  -min_distance to 0 
number_of_extra_samples = min_distance/(days+min_distance)*runs

extra_damp_high =[]
for i in range(int(number_of_extra_samples)):
    damping_days=(random.sample(range(min_damping_day,days), number_of_dampings))
    damping_days.append(random.randint(days, days+min_distance))
    extra_damp_high.append(sorted(damping_days))


extra_damp_low=[]
for i in range(int(number_of_extra_samples)):
    damping_days=(random.sample(range(min_damping_day,days), number_of_dampings))
    damping_days.append(random.randint(min_damping_day-min_distance, min_damping_day))
    extra_damp_low.append(sorted(damping_days))

additional_array = extra_damp_high+extra_damp_low


#calculate difference
differences = []
indices = []
for days_array, index in zip(additional_array, np.arange(0, len(additional_array)) ): 
    differences.append(np.diff(days_array).tolist())
    indices.append(index)

number_of_differences = len(differences[0])

res = [[i, d] for [i,d] in zip(indices, differences)  if len([j for j in d if j >= min_distance]) == number_of_differences]

#res = [[i, d] for [i,d] in zip(indices, differences)  if len([j for j in d if j >= min_distance]) == number_of_differences]
df_diff = pd.DataFrame(data = res, columns = ['index', 'diff'])

cleaned_additional = []
for index in sorted(df_diff['index'], reverse=True):
    cleaned_additional.append(additional_array[index])


#indices = np.arange(0, len(cleaned_additional))

#cleaned_indx =  [[i, d] for [i,d] in zip(indices, cleaned_additional)  if len([j for j in d if (j >= min_damping_day & j <=days)]) == number_of_dampings]


final_additional = []
for x in cleaned_additional : 
    keep = [j for j in x if j> min_damping_day and j <=days]
    if len(keep) == number_of_dampings:
        final_additional.append(keep)

# add additonals to originals 
final_array = cleaned_all_dampings+ final_additional

# only keep 10.000 of the list, to make it comparabel to the example before 
#final_array = final_array[:1000000]
final_array = random.sample(final_array, 1000000)



#plots for random sampling 
name = '_random'
plot_dist(all_dampings,days, name)
plot_differences(all_dampings, name)
plot_boxplot(all_dampings, name)
plot_maxminmean(all_dampings, name)

# plots for random sampling and then deleting 
name = '_deleted_new'
plot_dist(final_array, days, name)
plot_differences(final_array, name)
plot_boxplot(final_array, name)
plot_maxminmean(final_array, name)
