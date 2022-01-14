import os
from _thread import start_new_thread


n_runs = 40
samples = 300#000
n_dampings = 3



global marker
marker = 0

length = str(int(samples / n_runs))


def run(i):
    print(i)
    os.system("python3 ./multiple_data_generator.py " + length + " " + str(i) + " " + str(n_dampings))
    global marker
    marker += 1


for i in range(n_runs):
    start_new_thread(run,(i,))

while marker != n_runs:
    pass

csvfiles =["data/epi_data_"+str(i)+".csv" for i in range(n_runs)]
files_string =""
for item in csvfiles:
    files_string += item+" "
os.system("cat " + files_string + "> data/test.csv")
os.system("rm "+files_string)
