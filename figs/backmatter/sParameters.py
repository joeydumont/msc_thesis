# ------------------- Information --------------------- #
# Author:	Joey Dumont <joey.dumont@gmail.com>	#
# Date created:	October 7th, 2013			#
# Date mod.	October 7th, 2013			#
# Description:	We plot the scattering parameters of	#
#		RF-21 fibre design. Experimental and	#
#		theoretical.				#
# ----------------------------------------------------- #

# --------------- Modules Importation ----------------- #
from pylab import *
from matplotlib.ticker import AutoMinorLocator

# Setting the rc parameters.
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage[charter]{mathdesign}'
rcParams['font.size'] = 10
rcParams['legend.numpoints'] = 3

# ----------------- Data Importation ------------------ #
S11 = loadtxt("RF21-S11.dat")
Ssim = loadtxt("RF21-sParameters-Ag4.csv",skiprows=1,delimiter=",")
S22 = loadtxt("RF21-S22.dat")
S12 = loadtxt("RF21-S21.txt", skiprows=2)

# ---------------- Data Manipulation ------------------ #
# -- We plot the data as a function of normalized values
c=3.0e8/10e9
sz=32.0e-3
dz=28.0e-3
P=sz+dz
dl=27.0e-3
dr=55.0e-3


# -- We use two axes to plot data.
freqmin = min(S11[0,0],Ssim[0,0],S12[0,0],S22[0,0])
freqmax = max(S11[-1,0],Ssim[-1,0],S12[-1,0],S22[-1,0])
Smin = min(min(S11[:,1]),min(Ssim[:,1]),min(S12[:,1]),min(S22[:,1]))
Smax = max(max(S11[:,1]),max(Ssim[:,1]),max(S12[:,1]),max(S22[:,1]))

# ------------------ Plotting data -------------------- #
fig1 = figure(figsize=(7,3))
ax1 = fig1.add_subplot(111)

ax1.plot(Ssim[:,0],Ssim[:,1],"b--") 
ax1.plot(Ssim[:,0],Ssim[:,2],"k--")


minorLocator   = AutoMinorLocator()
ax1.xaxis.set_minor_locator(minorLocator)

ax1.set_xlabel("Frequency (GHz)")
ax1.set_ylabel("$S$-parameters (dB)")
xlim((0.1,5))

fig1.savefig("sParametersRF21sim.pdf", bbox_inches='tight')


ax1.plot(S11[:,0],S11[:,1], "b-", label="$S_{11}$")
ax1.plot(S22[:,0],S22[:,1], "r-", label="$S_{22}$")
ax1.plot(S12[:,0],S12[:,1], "k-", label="$S_{12}$")
ax1.legend(loc=0)

ax2 = ax1.twiny()
ax2.set_xlabel(r"$\lambda_f/P$")
ax2.set_xticks(ax1.get_xticks())
def tick_function(X):
  V = c/(P*X)
  return ["%.3f" % z for z in V]

ax2.set_xticklabels(tick_function(ax1.get_xticks()))
#show()
fig1.savefig("sParametersRF21.pdf", bbox_inches="tight")
