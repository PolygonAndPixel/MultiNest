####################evaluated data - only lowest likelihood per iteration#######
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib
import sys

cm = plt.cm.get_cmap('jet')

print "Using data from " + sys.argv[1]
print "Save picture as " + sys.argv[2]
text = ""
for i in range(3, len(sys.argv)):
    text = text + " " + sys.argv[i]
print "Title of graph is " + text

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

with open(sys.argv[1], 'r') as f:
    x = []
    y = []
    likelihood = []
    cluster = []
    for line in f:
        data = line.split()
        x.append(float(data[0]))
        y.append(float(data[1]))
        likelihood.append(float(data[2]))
        cluster = int(data[4])

# For a better color with Himmelblau's function, use: vmin=-100
sc = ax.scatter(x, y, likelihood, 'r', c=likelihood, cmap=cm, s=5, edgecolors='none')
ax.set_xlabel('X')
ax.set_ylabel('Y')
plt.suptitle(text)
fig.colorbar(sc)
sc.colorbar.set_label('Evaluated data')
plt.savefig(sys.argv[2], dpi=300)