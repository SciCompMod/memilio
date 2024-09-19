from matplotlib.ticker import ScalarFormatter, LogLocator
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
fontsize = 24

# First plot: Runtime vs. Number of Processors on
# processors = [1, 2, 4, 8]  # Number of processors
# runtime = [[0.2, 0.3, 0.6, 1.1],  # 0.025 mio agents per processor
#                 [0.6, 0.8, 1.2, 2.3],  # 0.05 mio agents per processor
#                 [1.3, 1.7, 2.4, 4.9],  # 0.1 mio agents per processor
#                 [2.7, 3.5, 5.3, 11]]  # 0.2 mio agents per processor
processors = [1, 2, 4, 8, 16, 32]  # Number of processors
runtime = [ [0.9, 1.3, 1.5, 1.9, 2.6, 3.9],  # 0.025 mio agents per processor
            [2.5, 2.8, 3.0, 3.9, 5.4, 7.0],  # 0.05 mio agents per processor
            [5.2, 5.8, 6.3, 7.9, 10.4, 14.7],  # 0.1 mio agents per processor
            [10.6, 11.1, 12.8, 15.6, 21.1, 32.1]]   # 0.2 mio agents per processor
runtime = np.array(runtime)  # Convert to seconds
fig = plt.figure(figsize=(12, 9))
for i in range(len(runtime)):
    plt.plot(processors, runtime[i], marker='o', linewidth=3)

plt.legend(['25k agents per processor', '50k agents per processor', '100k agents per processor', '200k agents per processor'], fontsize=fontsize-4)
plt.xlabel('Number of Processors', fontsize=fontsize)
plt.ylabel('Runtime (seconds)', fontsize=fontsize)
plt.title('Node level weak scaling', fontsize=fontsize+4)
plt.yscale('log')
plt.tick_params(axis='both', which='major', labelsize=fontsize-4)
plt.tick_params(axis='both', which='minor', labelsize=fontsize-4)
plt.grid(True)
plt.show()

# Second plot: Speedup vs. Number of Processors
# time for 6400k agents on x processors
# processors = [1, 2, 4, 8, 16, 32]  # Number of processors
# time_800k = [95, 49, 32, 19, 12.1, 8.3]  # Time in seconds
# time_6400k = [741, 407, 220, 143, 93, 59]  # Time in seconds
# speedup = [time_6400k[0] / t for t in time_6400k]
# ideal_speedup = [p for p in processors]

# plt.figure()
# plt.plot(processors, speedup, marker='o')
# plt.plot(processors, ideal_speedup, marker='o')
# plt.xlabel('Number of Processors')
# plt.ylabel('Speedup')
# plt.title('Speedup vs. Number of Processors for 6.4mio Agents')
# plt.legend(['Speedup', 'Ideal Speedup'])
# plt.grid(True)
# plt.show()



# Third plot: CPU Time vs. Population Size on 8 processors
population_size = [25, 50, 100, 200, 400, 800, 1600, 3200, 6400]  # Population size in thousand
population_size = [p * 1000 for p in population_size]
cpu_time_real = np.array([0.14, 0.36, 0.9, 1.8, 3.9, 7.92, 15.61, 31.2, 65.5])*(1/120.0)  # CPU time in seconds
cpu_time_real_1core = np.array([0.9, 2.51, 5.1, 10.4, 21.1, 42.5, 85.2, 173.2, 349.2])*(1/120.0)  # CPU time in seconds

# we also need the average number of seconds per 1000 agents
cpu_time_per_1000 = [t / p for t, p in zip(cpu_time_real, population_size)]
average_cpu_time_per_1000 = sum(cpu_time_per_1000) / len(cpu_time_per_1000)
# calculate the CPU time with the average for above population sizes
cpu_time_av = [p * average_cpu_time_per_1000 for p in population_size]

# Fourth plot: Memory Usage vs. Population Size
population_size_memory = [25, 50, 100, 200, 400, 800, 1600, 3200, 6400]  # Population size in thousand
population_size_memory = [p * 1000 for p in population_size_memory]
memory_usage = [10, 22, 50, 111, 228, 455, 984, 1910, 3780]  # Memory usage in MB
# we also need the average memory usage per 1000 agents
average_memory_time_per_1000 = [m / p for m, p in zip(memory_usage, population_size_memory)]
average_memory_usage = sum(average_memory_time_per_1000) / len(average_memory_time_per_1000)
mem_av = [p * average_memory_usage for p in population_size_memory]


# Plot for CPU Time vs. Population Size
plt.figure(figsize=(12, 9))
plt.plot(population_size, cpu_time_real, marker='o', linewidth=3)
plt.plot(population_size, cpu_time_av, linewidth=3, linestyle='dashed', color='black', label='Average CPU time')
plt.yscale('log')
plt.xscale('log')
plt.legend(['Runtime Simulation', 'Linear Scaling: ' +'%.2E' % Decimal(average_cpu_time_per_1000)+' Seconds / Agent'], fontsize=fontsize)
plt.xlabel('Population Size', fontsize=fontsize)
plt.ylabel('Runtime (seconds)', fontsize=fontsize)
plt.tick_params(axis='both', which='major', labelsize=fontsize-4)
plt.gca().yaxis.set_major_formatter(ScalarFormatter())
plt.gca().yaxis.set_minor_formatter(ScalarFormatter())
plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, numticks=15))
plt.xticks([10000, 100000, 1000000, 10000000])
plt.gca().get_xaxis().set_major_locator(LogLocator(base=10.0, numticks=15))
plt.gca().get_xaxis().set_minor_locator(LogLocator(base=10.0, numticks=15))
plt.title('Runtime Scaling', fontsize=fontsize+4)
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot for Memory Usage vs. Population Size
plt.figure(figsize=(12, 9))
plt.plot(population_size_memory, memory_usage, marker='o', linewidth=3)
plt.plot(population_size_memory, mem_av, linewidth=3, linestyle='dashed', color='black', label='Average Memory Usage')
plt.legend(['Memory Usage Simulation', 'Linear Scaling: ' + '%.2E' % Decimal(average_memory_usage)+' MB / Agent'], fontsize=fontsize)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Population Size', fontsize=fontsize)
plt.ylabel('Memory Usage (MB)', fontsize=fontsize)
plt.title('Memory Usage Scaling', fontsize=fontsize+4)
plt.tick_params(axis='both', which='major', labelsize=fontsize-4)
plt.gca().yaxis.set_major_formatter(ScalarFormatter())
plt.gca().yaxis.set_minor_formatter(ScalarFormatter())
plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, numticks=15))
plt.xticks([10000, 100000, 1000000, 10000000])
plt.gca().get_xaxis().set_major_locator(LogLocator(base=10.0, numticks=15))
plt.gca().get_xaxis().set_minor_locator(LogLocator(base=10.0, numticks=15))
plt.grid(True)
plt.tight_layout()
plt.show()
