import matplotlib.pyplot as plt

# First plot: Runtime vs. Number of Processors on
# processors = [1, 2, 4, 8]  # Number of processors
# runtime = [[0.6, 0.8, 1.1, 2.2],  # 0.025 mio agents per processor
#                 [1.4, 1.8, 2.5, 4.5],  # 0.05 mio agents per processor
#                 [3.3, 4.8, 6.4, 9.5]]  # 0.1 mio agents per processor
processors = [1, 2, 4, 8, 16, 32]  # Number of processors
runtime = [ [02.0, 02.9, 04.5, 05.3, 06.1, 08.2],  # 0.025 mio agents per processor
            [05.7, 06.2, 08.9, 10.0, 12.1, 16.7],  # 0.05 mio agents per processor
            [11.2, 12.5, 16.8, 19.0, 24.3, 32.1],  # 0.1 mio agents per processor
            [23.6, 24.9, 32.1, 37.5, 47.9, 63.4]]   # 0.2 mio agents per processor

plt.figure()
for i in range(len(runtime)):
    plt.plot(processors, runtime[i], marker='o')
plt.legend(['25k agents per processor', '50k agents per processor', '100k agents per processor', '200k agents per processor'])
plt.xlabel('Number of Processors')
plt.ylabel('Runtime (seconds)')
plt.title('Runtime vs. Number of Processors for 25/50/100/200 k Agents per Processor for 10 days')

plt.grid(True)
plt.show()

# Second plot: Speedup vs. Number of Processors
# time for 6400k agents on x processors
processors = [1, 2, 4, 8, 16, 32]  # Number of processors
time_800k = [95, 49, 32, 19, 12.1, 8.3]  # Time in seconds
time_6400k = [741, 407, 220, 143, 93, 59]  # Time in seconds
speedup = [time_6400k[0] / t for t in time_6400k]
ideal_speedup = [p for p in processors]

plt.figure()
plt.plot(processors, speedup, marker='o')
plt.plot(processors, ideal_speedup, marker='o')
plt.xlabel('Number of Processors')
plt.ylabel('Speedup')
plt.title('Speedup vs. Number of Processors for 6.4mio Agents')
plt.legend(['Speedup', 'Ideal Speedup'])
plt.grid(True)
plt.show()

# Third plot: CPU Time vs. Population Size on 8 processors
population_size = [25, 50, 100, 200, 400, 800, 1600, 3200,6400]  # Population size in thousand
population_size = [p * 1000 for p in population_size]
cpu_time_real = [0.32,0.84,2.5,5.3,10.1,19.2,37.5,72.9,148]  # CPU time in seconds
# we also need the avarage number of seconds per 1000 agents
cpu_time_per_1000 = [t / p for t, p in zip(cpu_time_real, population_size)]
avarage_cpu_time_per_1000 = sum(cpu_time_per_1000) / len(cpu_time_per_1000)
#calculate the CP time with the avarage for above population sizes
cpu_time_av = [p * avarage_cpu_time_per_1000 for p in population_size]


plt.figure()
#logarithmic scale
plt.yscale('log')
plt.xscale('log')
plt.plot(population_size, cpu_time_real, marker='o')
plt.plot(population_size, cpu_time_av, marker='x', alpha=0.5)
plt.legend(['Real CPU Time', 'Constant CPU time = '+str(avarage_cpu_time_per_1000)+' * Population Size'])
plt.xlabel('Population Size')
plt.ylabel('CPU Time (seconds)')
plt.title('CPU Time vs. Population Size in thousand')
plt.grid(True)
plt.show()

# Fourth plot: Memory Usage vs. Population Size
population_size = [25, 50, 100, 200, 400, 800, 1600, 3200]  # Population size in thousand
population_size = [p * 1000 for p in population_size]
memory_usage = [10,22, 50, 111, 228,455,900,1777]  # Memory usage in MB 

plt.figure()
plt.plot(population_size, memory_usage, marker='o')
plt.xlabel('Population Size')
plt.ylabel('Memory Usage (MB)')
plt.title('Memory Usage vs. Population Size')
plt.grid(True)
plt.show()