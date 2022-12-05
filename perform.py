import matplotlib.pyplot as plt
import numpy as np

barWidth = 0.25

input_size = [ 100000, 500000, 1000000 ]
serial = [ 0.099610, 0.553899, 14.190338 ]
parallel = [ 0.147379, 0.914620, 4.430356 ]

br1 = np.arange(len(input_size))
br2 = [ x + barWidth for x in br1 ]

plt.bar(br1, serial, color='r', width=barWidth, edgecolor='grey', label='Serial')
plt.bar(br2, parallel, color='g', width=barWidth, edgecolor='grey', label='Parallel')

plt.title('Execution Time by Input Size (nprocs=4)', fontweight='bold', fontsize=17)
plt.xlabel('Input Size (integer count)', fontweight='bold', fontsize=15)
plt.ylabel('Time Elapsed (seconds)', fontweight='bold', fontsize=15)
plt.xticks([r + barWidth for r in range(len(input_size))], input_size)
plt.legend()

plt.savefig('results.out.png')
