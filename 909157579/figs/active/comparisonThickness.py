# -------------------- Information -------------------- #
# Author:       Joey Dumont <joey.dumont@gmail.com>     #
# Date created: Feb. 13th, 2014							#
# Date mod.:    Feb. 13th, 2014							#
# Description:  We plot the thickness of the metallic	#
#				layers for different values of the 		#
#				concentration via the Bruggeman model. 	#
#				We compare the results for the bulk		#
#				silver conductivity and the thickness	#
#				corrected one (Fuchs-Sondheimer). 		#
# ----------------------------------------------------- #

# --------------- Modules Importation ----------------- #
# Importing numerical analysis package and
# morgenstemning color map.
from pylab import *
from numpy import roots
from matplotlib import cm
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
rcParams['legend.numpoints'] = 5

# -------------- Constants Declaration ---------------- #
sigma1 = 6.3e7
sigma2 = 6.79e6
a = 99.646e-6
tInner = 0.354e-6
L = 24.372e-2
Rinner = 80.5 
pi = 3.14159

lamb = 40e-9

# --------------- Function Definition ----------------- #

def effectiveConductivity(s1,s2,d,delta):
	q = (d*delta-1)*s1+(d*(1-delta)-1)*s2
	return 1.0/(2*d-2)*(q+sqrt(q*q+4*(d-1)*s1*s2))

def thicknessBulkConductivity(radShell,lengthFibre,BulkCond,Resistance):
	return -radShell+sqrt(radShell*radShell+lengthFibre/(pi*BulkCond*Resistance))

def thicknessFuchsSondheimer(radShell, lengthFibre, BulkCond, Resistance, Specularity, Ell):
	p = [1, 2*radShell, -lengthFibre/(pi*Resistance*BulkCond), -3*lengthFibre*Ell*(1-Specularity)/(pi*Resistance*BulkCond)]
	return roots(p)[1]

def conductivityFuchsSondheimer(BulkCond,thick, ell,Specularity):
	return BulkCond/(1+3*ell/(8*thick)*(1-Specularity))

# ------------------ MAIN FUNCTION -------------------- #
#-- We test the trigonometric for the cubic function. 
aR = 1
bR = 2*a
cR = -L/(sigma1*pi*Rinner)
dR = -3.0*L*lamb/(8.0*sigma1*Rinner)

pR = (3.0*aR*cR-bR*bR)/(3.0*aR)
qR = (2.0*bR**3-9.0*aR*bR*cR+27.0*aR*aR*dR)/(27.0*aR**3)

k=0
t0 = 2.0*sqrt(-pR/3.0)*cos(arccos(3.0*qR/(2.0*pR)*sqrt(-3.0/pR))/3.0-k*2.0*pi/3.0)
t = t0-bR/(3.0*aR)
print(t)
print(sigma1/(1.0+3.0*lamb/(8.0*t)))

#-- We plot the differences.
size = 100
concentration = linspace(1.0,0.0,size)


fig1 = figure(figsize=(6,3))
ax1 = fig1.add_subplot(111)
results = zeros((size,7))
for i in range(size):
	# Compute effective cond.
	results[i,0] = concentration[i]
	results[i,1] = effectiveConductivity(sigma1,sigma2,3,results[i,0])
	results[i,2] = thicknessBulkConductivity(a,L,results[i,1],Rinner)
	results[i,3] = thicknessFuchsSondheimer(a,L,results[i,1],Rinner,0,lamb)
	results[i,4] = conductivityFuchsSondheimer(results[i,1],results[i,3],lamb,0)

print(results[0,2])
	

ax1.set_xlabel(r"Silver concentration ($\delta_1$)")




ax2 = ax1.twinx()
ax2.set_ylabel(r"$\sigma/\sigma_0$")
ax2.plot(results[:,0], results[:,1]/sigma1, 'k-', label=r"Bruggeman")
ax2.plot(results[:,0], results[:,4]/sigma1, 'k--', label=r"Bruggeman+FS")
ax1.plot(results[:,0], results[:,2]/results[0,2], 'b-', label=r"Bruggeman")
ax1.plot(results[:,0], results[:,3]/results[0,2], 'b--', label=r"Bruggeman+FS")
ax1.set_ylabel("$t/t_0$", color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')

ax1.legend(loc=2)

savefig("comparisonThickness.pdf", bbox_inches='tight')

# ------------- Additional Computations --------------- #
def thicknessFuchsSondheimerInner(radShell, lengthFibre, BulkCond, Resistance, Specularity, Ell):
	p = [1, -2*radShell, lengthFibre/(pi*Resistance*BulkCond), 3*lengthFibre*Ell*(1-Specularity)/(4*pi*Resistance*BulkCond)]
	return roots(p)

def thicknessFuchsSondheimerOuter(radShell, lengthFibre, BulkCond, Resistance, Specularity, Ell):
	p = [1, 2*radShell, -lengthFibre/(pi*Resistance*BulkCond), -3*lengthFibre*Ell*(1-Specularity)/(4*pi*Resistance*BulkCond)]
	return roots(p)

innerThickness = thicknessFuchsSondheimerInner(100e-6, 24.372e-2,sigma1, 80.5, 0, lamb)
innerConductiv = sigma1/(1+3*lamb/(8*innerThickness))
print(innerThickness)
print(innerConductiv)

outerThickness = thicknessFuchsSondheimerOuter(397e-6, 24.372e-2, sigma1, 59.8, 0, lamb)
outerConductiv = sigma1/(1+3*lamb/(8*outerThickness))
print(outerThickness)
print(outerConductiv)