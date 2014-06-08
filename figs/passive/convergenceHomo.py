# -------------------- Information -------------------- #
# Author:       Joey Dumont <joey.dumont@gmail.com>     #
# Date created: Mar. 22nd, 2014							#
# Date mod.:    Jun. 6th, 2014							#
# Description:  We show the convergence properties of	#
#				SQA in the case of an homogeneous disk. #
# ----------------------------------------------------- #

# --------------- Modules Importation ----------------- #
# Importing numerical analysis package.
from pylab import *
from numpy import polyfit

# Setting the rc parameters.
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage[charter]{mathdesign}'
rcParams['font.size'] = 10
#rcParams['axes.labelsize'] = '10'
#rcParams['xtick.labelsize'] = '10'
#rcParams['ytick.labelsize'] = '10'
rcParams['legend.numpoints'] = 3

# ------------------ Loading Data --------------------- #
X = loadtxt("convergenceHomo.dat")

# -------------------- Plot Data ---------------------- #
fig1 = figure(figsize=(7,3))
k = array([2.5, 5.0, 10.0, 20.0, 40.0, 80.0, 120.0])

for i in range(int(X.shape[1]/2)):
  plot(X[:,2*i], X[:,2*i+1], label="$kR_0=%2.2f$" %(k[i]))

  # We elide part of the data.
  imin = int(X.shape[0]/4)
  print(polyfit(log(X[imin:-imin,2*i]),log(X[imin:-imin,2*i+1]), 1))

xlabel("$kR_02\epsilon$")
ylabel("Maximum error")
gca().invert_xaxis()
gca().set_yscale('log')
gca().set_xscale('log')
grid(True,alpha=0.5)
xlim((10,1e-4))
#ylim((1e-15,1e-3))
legend(loc=0, prop={'size':7})

savefig("convergenceHomo.pdf", bbox_inches='tight')