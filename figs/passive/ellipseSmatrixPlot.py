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
import morgenstemning as mrg
ms, msi = mrg.morgenstemning()

# Setting the rc parameters.
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage[charter]{mathdesign}'
rcParams['font.size'] = 10
rcParams['axes.labelsize'] = '10'
rcParams['xtick.labelsize'] = '10'
rcParams['ytick.labelsize'] = '10'
rcParams['legend.numpoints'] = 3

# ------------------ Loading Data --------------------- #
x = np.loadtxt("sMatrixEllipseMag.dat")
y = np.arange(-33,34,1)

# -------------------- Plot Data ---------------------- #
fig = plt.figure(figsize=(3,3))
ax = fig.add_subplot(111)
ax.set_xlim((-33,33))
ax.set_ylim((-33,33))
ax.invert_yaxis()
plt.pcolor(y,y,x, cmap=msi)
plt.savefig("sMatrixEllipseMag.pdf", bbox_inches='tight')