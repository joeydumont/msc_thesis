# ------------------- Information --------------------- #
# Author:	Joey Dumont <joey.dumont@gmail.com>	#
# Date created:	October 18th, 2013			#
# Date mod.	October 18th, 2013			#
# Description:	We plot the times it took to compute	#
#		sets of Wigner symbols of different	#
#		sizes.					#
# ----------------------------------------------------- #

# --------------- Modules Importation ----------------- #
from pylab import *
from matplotlib.ticker import AutoMinorLocator

# ----------------- Data Importation ------------------ #
prec = loadtxt("precisionSph.dat")

# ------------------ Plotting data -------------------- #
fig1 = figure(figsize=(7,3))
ax1 = fig1.add_subplot(111)

ax1.plot(prec[:,0],prec[:,1], 'b-')
ax1.plot(prec[:,0],prec[:,2], 'r')
ax1.plot(prec[:,0],prec[:,3], 'k')

minorLocator   = AutoMinorLocator()
ax1.xaxis.set_minor_locator(minorLocator)

ax1.set_xlabel(r"$\ell$")
ax1.set_ylabel("Error")
ax1.set_yscale('log')

fig1.savefig("SphPrecision.pdf", bbox_inches="tight")