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
times = loadtxt("times.dat")
timesF = loadtxt("timesF.dat")

prec = loadtxt("precision.dat")
precF = loadtxt("precisionF.dat")

# ---------------- Data Manipulation ------------------ #
# -- We plot the data in milliseconds.
times[:,1]*=1000
timesF[:,1]*=1000

# -- We program our histrogram. 
def histogram(X, nbins):
  dx = (max(X[:,0])-min(X[:,0]))/(nbins)
  centerBins = [min(X[:,0])+(n-0.5)*dx for n in range(1,nbins+1)]
  avgs = zeros((nbins))
  for i in range(nbins):
    Y = X[X[:,0]<(centerBins[i]+dx/2),:]
    Z = Y[Y[:,0]>(centerBins[i]-dx/2),1]
    avgs[i] = mean(Z)
  return centerBins, avgs, dx

# ------------------ Plotting data -------------------- #
fig1 = figure(figsize=(3,3))
ax1 = fig1.add_subplot(111)

centerBins, avgs, dx = histogram(times,10)
centerBinsF,avgsF,dxF=histogram(timesF,10)

ax1.bar(centerBins,avgs, width=dx/2)
ax1.bar(centerBinsF,avgsF,width=dxF/2, color='r')
ax1.plot(times[:,0],times[:,1], 'b-', alpha=0.75, label="C++")
ax1.plot(timesF[:,0],timesF[:,1], 'r', alpha=0.75, label="Fortran")
legend(loc=0,fancybox=True,shadow=True)

xlim(times[0,0], times[-1,0])
#ylim(0,150)

minorLocator   = AutoMinorLocator()
ax1.xaxis.set_minor_locator(minorLocator)

ax1.set_xlabel("Size")
ax1.set_ylabel("Time (ms)")

#show()
fig1.savefig("wignerTimes.pdf", bbox_inches="tight")

fig2 = figure(figsize=(3,3))
ax2 = fig2.add_subplot(111)

ax2.plot(prec[:,0],prec[:,1],'b-')
ax2.plot(precF[:,0],precF[:,1], 'r-')
ax2.set_yscale('log')

ax1.set_xlabel("$m_3$")
ax2.set_ylabel("Error")

show()
fig2.savefig("wignerPrecision.pdf", bbox_inches='tight')