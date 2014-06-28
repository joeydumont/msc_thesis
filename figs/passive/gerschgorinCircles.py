# -------------------- Information -------------------- #
# Author:       Joey Dumont <joey.dumont@gmail.com>     #
# Date created: May 26th, 2014							#
# Date mod.:    May 27th, 2014							#
# Description:  We compute the Gerschgorin circles of	#
#				any given matrix, and then one in 		#
#				particular. 							#
# ----------------------------------------------------- #

# --------------- Modules Importation ----------------- #
# Importing numerical analysis package.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Setting the rc parameters.
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage[charter]{mathdesign}'
rcParams['font.size'] = 10
rcParams['axes.labelsize'] = '10'
rcParams['xtick.labelsize'] = '10'
rcParams['ytick.labelsize'] = '10'
rcParams['legend.numpoints'] = 3

def gerschgorinCircles(mat, mode="full"):
	# -- Check if input array is a square matrix.
	if mat.ndim != 2:
		print("Input must a 2D array, i.e. a matrix.")
		return 0

	if mat.shape[0]!=mat.shape[1]:
		print("Input must be a square matrix!")
		return 0

	# -- Compute the positions of the circles and the
	# -- associated row and column radii.
	centerPoints = np.zeros((0,2))
	radii = np.zeros((mat.shape[0],2))
	for i in range(mat.shape[0]):
		centerPoints = np.insert(centerPoints, len(centerPoints), [np.real(mat[i,i]), np.imag(mat[i,i])], axis=0)

		for j in range(mat.shape[0]):
			if j!=i:
				# Row radius
				radii[i,0] += np.abs(mat[i,j])

				# Column radius
				radii[i,1] += np.abs(mat[j,i])

	# -- Depending on the mode, select the proper slice of the data. 
	if mode=="row":
		radii = radii[:,0]
	elif mode=="col":
		radii = radii[:,1]
	elif mode=="full":
		radii = np.min(radii,axis=1)

	# -- Plot the circles.
	fig = plt.figure(figsize=(5,3))
	ax = fig.add_subplot(111)
	ax.set_xlim((-8,8))
	ax.set_ylim((-8*3/5,8*3/5))
	ax.spines['left'].set_position('zero')
	ax.spines['right'].set_color('none')
	ax.spines['bottom'].set_position('zero')
	ax.spines['top'].set_color('none')
	#ax.spines['left'].set_smart_bounds(True)
	#ax.spines['bottom'].set_smart_bounds(True)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks(np.arange(-8,9,1))
	ax.yaxis.set_ticks(np.arange(-4,5,1))
	print(centerPoints)
	print(radii)

	for i in range(mat.shape[0]):
		circ = plt.Circle((centerPoints[i,0], centerPoints[i,1]), radii[i], facecolor=(0.7,0.7,0.7,0.4), edgecolor="black")
		ax.add_patch(circ)


	return fig

A = np.array([[5, 1, 1], [0, 6, 1], [1, 0, -5]])
B = np.array([[1,1/2,1/3], [1/2,1/3,1/4],[1/3,1/4,1/5]])
fig = gerschgorinCircles(A, "row")
plt.savefig("gerschgorin-row.pdf", bbox_inches='tight')

fig = gerschgorinCircles(A, "col")
plt.savefig("gerschgorin-col.pdf", bbox_inches='tight')

fig = gerschgorinCircles(A, "full")
plt.savefig("gerschgorin-full.pdf", bbox_inches='tight')