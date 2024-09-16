from matplotlib.ticker import ScalarFormatter, LogLocator

import matplotlib.pyplot as plt
from decimal import Decimal
fontsize = 18

# First plot: Runtime vs. Number of Processors on
# processors = [1, 2, 4, 8]  # Number of processors
# runtime = [[0.6, 0.8, 1.1, 2.2],  # 0.025 mio agents per processor
#                 [1.4, 1.8, 2.5, 4.5],  # 0.05 mio agents per processor
#                 [3.3, 4.8, 6.4, 9.5]]  # 0.1 mio agents per processor
# processors = [1, 2, 4, 8, 16, 32]  # Number of processors
# runtime = [ [02.0, 02.9, 04.5, 05.3, 06.1, 08.2],  # 0.025 mio agents per processor
#             [05.7, 06.2, 08.9, 10.0, 12.1, 16.7],  # 0.05 mio agents per processor
#             [11.2, 12.5, 16.8, 19.0, 24.3, 32.1],  # 0.1 mio agents per processor
#             [23.6, 24.9, 32.1, 37.5, 47.9, 63.4]]   # 0.2 mio agents per processor
# fig = plt.figure(figsize=(12, 9))
# for i in range(len(runtime)):
#     plt.plot(processors, runtime[i], marker='o', linewidth=3)

# plt.legend(['25k agents per processor', '50k agents per processor', '100k agents per processor', '200k agents per processor'], fontsize=fontsize)
# plt.xlabel('Number of Processors', fontsize=fontsize)
# plt.ylabel('Runtime (seconds)', fontsize=fontsize)
# plt.title('Node level weak scaling', fontsize=fontsize+4)
# plt.grid(True)
# plt.show()

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



# # Third plot: CPU Time vs. Population Size on 8 processors
population_size = [25, 50, 100, 200, 400, 800, 1600, 3200, 6400]  # Population size in thousand
population_size = [p * 1000 for p in population_size]
cpu_time_real = [0.32, 0.84, 2.5, 5.3, 10.1, 19.2, 37.5, 72.9, 148]  # CPU time in seconds
# we also need the average number of seconds per 1000 agents
cpu_time_per_1000 = [t / p for t, p in zip(cpu_time_real, population_size)]
average_cpu_time_per_1000 = sum(cpu_time_per_1000) / len(cpu_time_per_1000)
# calculate the CPU time with the average for above population sizes
cpu_time_av = [p * average_cpu_time_per_1000 for p in population_size]

# Fourth plot: Memory Usage vs. Population Size
population_size_memory = [25, 50, 100, 200, 400, 800, 1600, 3200]  # Population size in thousand
population_size_memory = [p * 1000 for p in population_size_memory]
memory_usage = [10, 22, 50, 111, 228, 455, 900, 1777]  # Memory usage in MB
average_memory_usage = sum(memory_usage) / len(memory_usage) * 1e-6
mem_av = [p * average_memory_usage for p in population_size_memory]


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9))

# Plot for CPU Time vs. Population Size

ax1.plot(population_size, cpu_time_real, marker='o', linewidth=3)
ax1.plot(population_size, cpu_time_av, linewidth=3, linestyle='dashed', color='black', label='Average CPU time')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.legend(['Runtime Simulation', 'Linear Scaling:' +'%.2E' % Decimal(average_cpu_time_per_1000)+' Seconds / Agent'], fontsize=fontsize)
ax1.set_xlabel('Population Size',fontsize=fontsize)
ax1.set_ylabel('Runtime (seconds)',fontsize=fontsize)
ax1.tick_params(axis='both', which='major', labelsize=fontsize-4)
ax1.yaxis.set_major_formatter(ScalarFormatter())
ax1.yaxis.set_minor_formatter(ScalarFormatter())
ax1.yaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
ax1.yaxis.set_minor_locator(LogLocator(base=10.0, numticks=15))
# x ticks should be from 10000 to 10000000
ax1.set_xticks([10000, 100000, 1000000, 10000000])
ax1.get_xaxis().set_major_locator(LogLocator(base=10.0, numticks=15))
ax1.get_xaxis().set_minor_locator(LogLocator(base=10.0, numticks=15))
ax1.set_title('Runtime Scaling',fontsize=fontsize+4)
ax1.grid(True)

# Plot for Memory Usage vs. Population Size
ax2.plot(population_size_memory, memory_usage, marker='o', linewidth=3)
ax2.plot(population_size_memory, mem_av, linewidth=3, linestyle='dashed', color='black', label='Average Memory Usage')
ax2.legend(['Memory Usage Simulation', 'Linear Scaling:' + '%.2E' % Decimal(average_memory_usage)+' MB / Agent'], fontsize=fontsize)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('Population Size',fontsize=fontsize)
ax2.set_ylabel('Memory Usage (MB)',fontsize=fontsize)
ax2.set_title('Memory Usage Scaling',fontsize=fontsize+4)
ax2.tick_params(axis='both', which='major', labelsize=fontsize-4)
ax2.yaxis.set_major_formatter(ScalarFormatter())
ax2.yaxis.set_minor_formatter(ScalarFormatter())
ax2.yaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
ax2.yaxis.set_minor_locator(LogLocator(base=10.0, numticks=15))
# x ticks should be from 10000 to 10000000
ax2.set_xticks([10000, 100000, 1000000, 10000000])
ax2.get_xaxis().set_major_locator(LogLocator(base=10.0, numticks=15))
ax2.get_xaxis().set_minor_locator(LogLocator(base=10.0, numticks=15))
ax2.grid(True)

plt.tight_layout()
plt.show()