'''Script that reads in a current file, defined as infile, and plots the values. Missing values (-99) should resemble Hawaiian Archipelago if correct. Output needs rotating 90 degrees.'''

# Read in current file
infile = "hawaii8km_/hawaii8km_09_05_01_100m_v.txt"
testdata = dict([ (y, [ [0.0]*(251 + 1) for x in xrange(438 + 1)], ) for y in xrange(0, 0 + 1) ])
with open(infile, "r") as file3:
    for j in xrange(1, 251+1):
	for i in xrange(1, 438+1):
            testdata[0][i][j] = float(file3.readline())

# Plot figure with missing values
import matplotlib.pylab as plt
fig = plt.figure()
ax=fig.add_subplot(1,1,1)
ax.set_aspect('equal')
plt.imshow(testdata[0], interpolation='nearest', cmap=plt.cm.ocean)
plt.colorbar()
plt.show()

