tasks = 16384
# tasks = 1280
nodes = [1, 2, 4, 8, 16, 32, 64, 128, 168]

parallel_task = 0
for i in nodes:
    parallel_task += tasks / i

print(round(parallel_task + 0.49))