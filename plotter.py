import matplotlib.pyplot as plt;

input()
input()
data = [int(i) for i in input().strip().split(' ')]
label = [str(i) for i in data]

fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1])
ax.bar(label, data)
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
plt.autoscale(enable=True, axis='both', tight=True)
plt.savefig('data.out.png')
