import os
from _thread import start_new_thread


# number of different threads used to calculate simulations
num_threads = 10
num_samples = 200#500#000
num_dampings = 3
num_sim_days = 10
num_time_steps_between_sim_days = 1

final_dataset_path = "../data/county_epi_data_test.csv"
tmp_datasets_path = "../data/tmp_county_epi_data_"


global marker
marker = 0

num_sample_per_thread = str(int(num_samples / num_threads))


def run(i):
    print(i)
    os.system("python3 ./county_data_generator.py " + num_sample_per_thread + " " + tmp_datasets_path + str(i) + ".csv"  + " " + str(num_dampings) + " " + str(num_sim_days))
    global marker
    marker += 1


for i in range(num_threads):
    start_new_thread(run,(i,))

while marker != num_threads:
    pass

csvfiles =[tmp_datasets_path+str(i)+".csv" for i in range(num_threads)]
files_string =""
for item in csvfiles:
    files_string += item+" "
os.system("cat " + files_string + "> " + final_dataset_path)
os.system("rm "+files_string)
