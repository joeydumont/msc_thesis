# -------------------- Information -------------------- #
# Author:       Joey Dumont <joey.dumont@gmail.com>     #
# Date created: May 19th, 2014							#
# Date mod.:    May 19th, 2014							#
# Description:  We plot the scattering parameters of	#
#				fibre-antenna design as a function of	#
#				frequency. 								#
#				We then study the autocorrelation of	#
#				each signal and analyze their form:		#
#					- exponential;						#
#					- Lorentzian, etc. 					#
#				We also compute the Fourier transform	#
#				of each curve and see what it gives. 	#
# ----------------------------------------------------- #

# --------------- Modules Importation ----------------- #
# Importing numerical analysis package and
# morgenstemning color map.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams
from scipy import interpolate
from scipy.optimize import curve_fit
#from scipy.stats import pearsonr
#from numpy import convolve
import scipy.fftpack as fft
import morgenstemning as mrg
ms,msi = mrg.morgenstemning()

# Setting the rc parameters.
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = [r"\usepackage[charter]{mathdesign}"]
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 8
rcParams['axes.labelsize'] = 'large'
rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['legend.numpoints'] = 3
rcParams['figure.figsize'] = 7,3

# ---------------- Data Importation ------------------ #
# We import the data for each fibre and rescale the 
# frequencies to GHz.

# RF10
rf10S22 = np.loadtxt("../xpData/RF10-S22-OneEnd-2.txt", skiprows=1)
rf10S22[:,0] /= 1.0e9

# RF21
rf21S22 = np.loadtxt("../xpData/RF21-S22-OneEnd.txt", skiprows=1)
rf21S22[:,0] /= 1.0e9

# RF27
rf27 = np.loadtxt("../xpData/RF27-Sparam-01-5GHz.txt", skiprows=1)
rf27[:,0] /= 1.0e9

# RF29
rf29 = np.loadtxt("../xpData/RF29-S-param.txt", skiprows=1)
rf29[:,0] /= 1.0e9

# RF33
rf33 = np.loadtxt("../xpData/RF33-S-parameters.txt", skiprows=1)
rf33[:,0] /= 1.0e9

# ------------------ Plotting Data -------------------- #

# -- We plot the S-parameters of each fibre-antenna.

# RF10
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
plt.plot(rf10S22[:,0],rf10S22[:,1], 'r-', lw=1.5)
ax1.set_xlabel("Frequency (GHz)")
ax1.set_ylabel(r"$S_{22}$ (dB)")
ax1.set_ylim((-50,0))
ax1.grid(True)
ax1.text(0.01,0.93, "RF10", transform=ax1.transAxes, fontsize=10, bbox=dict(facecolor='gray',alpha=0.2,lw=0.0))

plt.savefig("RF10-sParameters.pdf", bbox_inches='tight')

# RF21
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
plt.plot(rf21S22[:,0], rf21S22[:,1], 'r-', lw=1.5)
ax2.set_xlabel("Frequency (GHz)")
ax2.set_ylabel(r"$S_{22}$ (dB)")
ax2.set_ylim((-50,0))
ax2.grid(True)
ax2.text(0.01,0.93, "RF21", transform=ax2.transAxes, fontsize=10, bbox=dict(facecolor='gray',alpha=0.2,lw=0.0))

plt.savefig("RF21-sParameters.pdf", bbox_inches='tight')

# RF27
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
plt.plot(rf27[:,0], rf27[:,1], 'b-', lw=1.5, label=r"$S_{11}$")
plt.plot(rf27[:,0], rf27[:,2], 'k-', lw=1.5, label=r"$S_{12}$")
plt.plot(rf27[:,0], rf27[:,4], 'r-', lw=1.5, label=r"$S_{22}$")
ax3.set_xlabel("Frequency (GHz)")
ax3.set_ylabel(r"$S$-parameters")
ax3.set_ylim((-50,0))
ax3.grid(True)
ax3.text(0.01,0.93, "RF27", transform=ax3.transAxes, fontsize=10, bbox=dict(facecolor='gray',alpha=0.2,lw=0.0))
ax3.legend(loc=0)

plt.savefig("RF27-sParameters.pdf", bbox_inches='tight')

# RF29
fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
plt.plot(rf29[:,0], rf29[:,1], 'b-', lw=1.5, label=r"$S_{11}$")
plt.plot(rf29[:,0], rf29[:,2], 'k-', lw=1.5, label=r"$S_{12}$")
plt.plot(rf29[:,0], rf29[:,4], 'r-', lw=1.5, label=r"$S_{22}$")
ax4.set_xlabel("Frequency (GHz)")
ax4.set_ylabel(r"$S$-parameters")
ax4.set_ylim((-50,0))
ax4.grid(True)
ax4.text(0.01,0.93, "RF29", transform=ax4.transAxes, fontsize=10, bbox=dict(facecolor='gray',alpha=0.2,lw=0.0))
ax4.legend(loc=0)

plt.savefig("RF29-sParameters.pdf", bbox_inches='tight')

# RF33
fig5 = plt.figure()
ax5 = fig5.add_subplot(111)
plt.plot(rf33[:,0], rf33[:,1], 'b-', lw=1.5, label=r"$S_{11}$")
plt.plot(rf33[:,0], rf33[:,2], 'k-', lw=1.5, label=r"$S_{12}$")
plt.plot(rf33[:,0], rf33[:,4], 'r-', lw=1.5, label=r"$S_{22}$")
ax5.set_xlabel("Frequency (GHz)")
ax5.set_ylabel(r"$S$-parameters")
ax5.set_ylim((-50,0))
ax5.grid(True)
ax5.text(0.01,0.93, "RF33", transform=ax5.transAxes, fontsize=10, bbox=dict(facecolor='gray',alpha=0.2,lw=0.0))
ax5.legend(loc=0)

plt.savefig("RF33-sParameters.pdf", bbox_inches='tight')

# ----------------- Data Manipulation ----------------- #
# We compute the autocorrelation of each data sets. 	#
# ----------------------------------------------------- #

# -- Function definition.
def lorentz(x,a,gamma):
	return a*0.5*gamma/(x*x+gamma**2/4.0)

def exponentiel(x, a):
	return exp(-a*x)

def gaussian(x,a):
	return exp(-a*x*x)

# -- We compute the autocorrelations.

# We compute the vector of delays.
delays = np.zeros((rf10S22.shape[0]))
for i in range(rf10S22.shape[0]):
	delays[i] = i*(rf10S22[1,0]-rf10S22[0,0])

halfDelay = len(delays)//2
pInit = [0.5, 0.5]

# RF10
meanRF10 = np.mean(rf10S22[:,1])
dataRF10 = rf10S22[:,1]-meanRF10
normRF10 = np.sum(dataRF10**2)
acorRF10 = np.correlate(dataRF10,dataRF10, "full")/normRF10
acorRF10 = acorRF10[len(acorRF10)/2:]
optLorentzRF10, covLorentzRF10 =  curve_fit(lorentz, delays[:halfDelay], acorRF10[:halfDelay], p0=pInit)

fig6 = plt.figure()
ax6 = fig6.add_subplot(111)
plt.plot(delays,acorRF10, 'r')
plt.plot(delays[:halfDelay],lorentz(delays[:halfDelay],optLorentzRF10[0],optLorentzRF10[1]), 
	color='r', linestyle='--', marker='^',markevery=4, label=r'$S_{22}$: $A=%.3g$, $\Gamma=%.3g$' %(optLorentzRF10[0], optLorentzRF10[1]))
ax6.set_xlabel("Frequency Shift $\Delta f$ (GHz)")
ax6.set_ylabel(r"Autocorrelation of $S_{22}$")
ax6.legend(loc=0)
ax6.text(0.01,0.93, "RF10", transform=ax6.transAxes, fontsize=10, bbox=dict(facecolor='gray',alpha=0.2,lw=0.0))

plt.savefig("RF10-autoCorrelation.pdf", bbox_inches='tight')

# RF21
delays = np.zeros((rf21S22.shape[0]))
for i in range(rf10S22.shape[0]):
	delays[i] = i*(rf21S22[1,0]-rf21S22[0,0])

halfDelay = len(delays)//2
pInit = [0.5, 0.5]

meanRF21 = np.mean(rf21S22[:,1])
dataRF21 = rf21S22[:,1]-meanRF21
normRF21 = np.sum(dataRF21**2)
acorRF21 = np.correlate(dataRF21,dataRF21, "full")/normRF21
acorRF21 = acorRF21[len(acorRF21)//2:]
optLorentzRF21, covLorentzRF21 = curve_fit(lorentz,delays[:halfDelay], acorRF21[:halfDelay], p0=pInit)

fig7 = plt.figure()
ax7 = fig7.add_subplot(111)
plt.plot(delays,acorRF21, 'r')
plt.plot(delays[:halfDelay],lorentz(delays[:halfDelay],optLorentzRF21[0],optLorentzRF21[1]), 
	color='r', linestyle='--', marker='^',markevery=4, label=r'$S_{22}$: $A=%.3g$, $\Gamma=%.3g$' %(optLorentzRF21[0], optLorentzRF21[1]))
ax7.set_xlabel("Frequency Shift $\Delta f$ (GHz)")
ax7.set_ylabel(r"Autocorrelation of $S_{22}$")
ax7.legend(loc=0)
ax7.text(0.01,0.93, "RF21", transform=ax7.transAxes, fontsize=10, bbox=dict(facecolor='gray',alpha=0.2,lw=0.0))

plt.savefig("RF21-autoCorrelation.pdf", bbox_inches='tight')

# RF27
delays = np.zeros((rf27.shape[0]))
for i in range(rf27.shape[0]):
	delays[i] = i*(rf27[1,0]-rf27[0,0])

halfDelay = len(delays)//2
autoCorrRF27 = np.zeros((len(delays),3))
curveParamRF27 = np.zeros((2,3))

j=0
for i in [1,2,4]:
	meanRF27 = np.mean(rf27[:,i])
	dataRF27 = rf27[:,i]-meanRF27
	normRF27 = np.sum(dataRF27**2)
	acorRF27 = np.correlate(dataRF27,dataRF27,"full")/normRF27
	acorRF27 = acorRF27[len(acorRF27)//2:]
	optLorentzRF27, covLorentzRF27 = curve_fit(lorentz, delays[:halfDelay], acorRF27[:halfDelay], p0=pInit)
	
	autoCorrRF27[:,j] = acorRF27
	curveParamRF27[:,j] = optLorentzRF27
	j += 1


fig8 = plt.figure()
ax8 = fig8.add_subplot(111)
plt.plot(delays,autoCorrRF27[:,0], 'b')
plt.plot(delays[:halfDelay],lorentz(delays[:halfDelay],curveParamRF27[0,0],curveParamRF27[1,0]), 
	color='b', linestyle='--', marker='o',markevery=8, label=r'$S_{11}$: $A=%.3g$, $\Gamma=%.3g$' %(curveParamRF27[0,0], curveParamRF27[1,0]))
plt.plot(delays,autoCorrRF27[:,1], 'k')
plt.plot(delays[:halfDelay],lorentz(delays[:halfDelay],curveParamRF27[0,1],curveParamRF27[1,1]), 
	color='k', linestyle='--', marker='d',markevery=8, label=r'$S_{12}$: $A=%.3g$, $\Gamma=%.3g$' %(curveParamRF27[0,1], curveParamRF27[1,1]))
plt.plot(delays,autoCorrRF27[:,2], 'r')
plt.plot(delays[:halfDelay],lorentz(delays[:halfDelay],curveParamRF27[0,2],curveParamRF27[1,2]), 
	color='r', linestyle='--', marker='^',markevery=8, label=r'$S_{22}$: $A=%.3g$, $\Gamma=%.3g$' %(curveParamRF27[0,2], curveParamRF27[1,2]))
ax8.set_xlabel("Frequency Shift $\Delta f$ (GHz)")
ax8.set_ylabel(r"Autocorrelation")
ax8.legend(loc=0)
ax8.text(0.01,0.93, "RF27", transform=ax8.transAxes, fontsize=10, bbox=dict(facecolor='gray',alpha=0.2,lw=0.0))

plt.savefig("RF27-autoCorrelation.pdf", bbox_inches='tight')

# RF29
delays = np.zeros((rf29.shape[0]))
for i in range(rf29.shape[0]):
	delays[i] = i*(rf29[1,0]-rf29[0,0])

halfDelay = len(delays)//2
autoCorrRF29 = np.zeros((len(delays),3))
curveParamRF29 = np.zeros((2,3))

j=0
for i in [1,2,4]:
	meanRF29 = np.mean(rf29[:,i])
	dataRF29 = rf29[:,i]-meanRF29
	normRF29 = np.sum(dataRF29**2)
	acorRF29 = np.correlate(dataRF29,dataRF29,"full")/normRF29
	acorRF29 = acorRF29[len(acorRF29)//2:]
	optLorentzRF29, covLorentzRF29 = curve_fit(lorentz, delays[:halfDelay], acorRF29[:halfDelay], p0=pInit)
	
	autoCorrRF29[:,j] = acorRF29
	curveParamRF29[:,j] = optLorentzRF29
	j += 1


fig9 = plt.figure()
ax9 = fig9.add_subplot(111)
plt.plot(delays,autoCorrRF29[:,0], 'b')
plt.plot(delays[:halfDelay],lorentz(delays[:halfDelay],curveParamRF29[0,0],curveParamRF29[1,0]), 
	color='b', linestyle='--', marker='o',markevery=8, label=r'$S_{11}$: $A=%.3g$, $\Gamma=%.3g$' %(curveParamRF29[0,0], curveParamRF29[1,0]))
plt.plot(delays,autoCorrRF29[:,1], 'k')
plt.plot(delays[:halfDelay],lorentz(delays[:halfDelay],curveParamRF29[0,1],curveParamRF29[1,1]), 
	color='k', linestyle='--', marker='d',markevery=8, label=r'$S_{12}$: $A=%.3g$, $\Gamma=%.3g$' %(curveParamRF29[0,1], curveParamRF29[1,1]))
plt.plot(delays,autoCorrRF29[:,2], 'r')
plt.plot(delays[:halfDelay],lorentz(delays[:halfDelay],curveParamRF29[0,2],curveParamRF29[1,2]), 
	color='r', linestyle='--', marker='^',markevery=8, label=r'$S_{22}$: $A=%.3g$, $\Gamma=%.3g$' %(curveParamRF29[0,2], curveParamRF29[1,2]))
ax9.set_xlabel("Frequency Shift $\Delta f$ (GHz)")
ax9.set_ylabel(r"Autocorrelation")
ax9.legend(loc=0)
ax9.text(0.01,0.93, "RF29", transform=ax9.transAxes, fontsize=10, bbox=dict(facecolor='gray',alpha=0.2,lw=0.0))

plt.savefig("RF29-autoCorrelation.pdf", bbox_inches='tight')

# RF33
delays = np.zeros((rf33.shape[0]))
for i in range(rf33.shape[0]):
	delays[i] = i*(rf33[1,0]-rf33[0,0])

halfDelay = len(delays)//2
autoCorrRF33 = np.zeros((len(delays),3))
curveParamRF33 = np.zeros((2,3))

j=0
for i in [1,2,4]:
	meanRF33 = np.mean(rf33[:,i])
	dataRF33 = rf33[:,i]-meanRF33
	normRF33 = np.sum(dataRF33**2)
	acorRF33 = np.correlate(dataRF33,dataRF33,"full")/normRF33
	acorRF33 = acorRF33[len(acorRF33)//2:]
	optLorentzRF33, covLorentzRF33 = curve_fit(lorentz, delays[:halfDelay], acorRF33[:halfDelay], p0=pInit)
	
	autoCorrRF33[:,j] = acorRF33
	curveParamRF33[:,j] = optLorentzRF33
	j += 1


fig10 = plt.figure()
ax10 = fig10.add_subplot(111)
plt.plot(delays,autoCorrRF33[:,0], 'b')
plt.plot(delays[:halfDelay],lorentz(delays[:halfDelay],curveParamRF33[0,0],curveParamRF33[1,0]), 
	color='b', linestyle='--', marker='o',markevery=8, label=r'$S_{11}$: $A=%.3g$, $\Gamma=%.3g$' %(curveParamRF33[0,0], curveParamRF33[1,0]))
plt.plot(delays,autoCorrRF33[:,1], 'k')
plt.plot(delays[:halfDelay],lorentz(delays[:halfDelay],curveParamRF33[0,1],curveParamRF33[1,1]), 
	color='k', linestyle='--', marker='d',markevery=8, label=r'$S_{12}$: $A=%.3g$, $\Gamma=%.3g$' %(curveParamRF33[0,1], curveParamRF33[1,1]))
plt.plot(delays,autoCorrRF33[:,2], 'r')
plt.plot(delays[:halfDelay],lorentz(delays[:halfDelay],curveParamRF33[0,2],curveParamRF33[1,2]), 
	color='r', linestyle='--', marker='^',markevery=8, label=r'$S_{22}$: $A=%.3g$, $\Gamma=%.3g$' %(curveParamRF33[0,2], curveParamRF33[1,2]))
ax10.set_xlabel("Frequency Shift $\Delta f$ (GHz)")
ax10.set_ylabel(r"Autocorrelation")
ax10.legend(loc=0)
ax10.text(0.01,0.93, "RF33", transform=ax10.transAxes, fontsize=10, bbox=dict(facecolor='gray',alpha=0.2,lw=0.0))

plt.savefig("RF33-autoCorrelation.pdf", bbox_inches='tight')


# ----------------- Data Manipulation ----------------- #
# We compute the Fourier transform of the data.			#
# ----------------------------------------------------- #
fig11 = plt.figure()
ax11 = fig11.add_subplot(111)
plotArgs = [dict(linestyle='-', color='b', label=r"$S_{11}$"),
			dict(linestyle='-', color='k', label=r"$S_{12}$"),
			dict(linestyle='-', color='r', label=r"$S_{22}$")]
j=0
for i in [1,2,4]:
	fftRF33 = fft.fftshift(fft.fft(rf33[:,i]))
	fftRF33 /= len(fftRF33)
	freqRF33 = fft.fftshift(fft.fftfreq(rf33.shape[0],d=rf33[1,0]-rf33[0,0]))
	plt.plot(freqRF33,abs(fftRF33), **plotArgs[j])
	j += 1

ax11.set_xlabel("FFT Frequency")
ax11.set_ylabel("Intensity of the FFT Spectrum")
ax11.set_xlim((0,20))
ax11.set_yscale('log')
ax11.legend(loc=0)
ax11.text(0.01,0.93, "RF33", transform=ax11.transAxes, fontsize=10, bbox=dict(facecolor='gray',alpha=0.2,lw=0.0))
ax11.set_xticks(np.arange(0,20.5,1.0))
ax11.grid(True)
plt.savefig("RF33-fft.pdf", bbox_inches='tight')