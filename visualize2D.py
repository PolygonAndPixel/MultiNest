import numpy as np
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


with open(sys.argv[1], 'r') as f:
    x = []
    y = []
    likelihood = []
    cluster = []
    for line in f:
        data = line.split()
        x.append(float(data[0]))
        y.append(float(data[1]))
        # Use for gaussian shell
        #likelihood.append(float(data[2]))
        # Use for Himmelblau and eggbox
        likelihood.append(float(data[2]))
# For a better color with Himmelblau's function, use: vmin=-100
sc = plt.scatter(x, y, c=likelihood, cmap=cm, edgecolors='none', s=5)   
plt.colorbar(sc)
plt.suptitle(text)
plt.xlabel('X')
plt.ylabel('Y')
llhLabel = 'Evaluated data' 
sc.colorbar.set_label(llhLabel)
plt.savefig(sys.argv[2], dpi=300)
#plt.show()
