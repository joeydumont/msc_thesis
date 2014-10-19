# ------------------- Information --------------------- #
# Author:		Joey Dumont <joey.dumont@gmail.com>		#
# Date created:	October 23rd, 2013						#
# Date mod.	:	October 23rd, 2013						#
# Description:	We compute the (complex) time delays	#
#				for the homoneous, circular cavity		#
#				as a function of the complex part of	#
#				potential. 								#
# ----------------------------------------------------- #

# --------------- Modules Importation ----------------- #
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigvals
from scipy.special import jn,jvp,hankel1,h1vp,hankel2,h2vp
from scipy.optimize import newton
from scipy.integrate import quadrature
import morgenstemning as mrg
mr, mri = mrg.morgenstemning()

# Setting the rc parameters.
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = [r"\usepackage[charter]{mathdesign}"]
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 10
plt.rcParams['legend.numpoints'] = 3

# --------------- Function Definition ----------------- #
def timeDelayEigenfunctionRe(r,m,nc,k):
	eps2 = np.imag(nc*nc)
	return np.real(-1j*k*eps2*2.0*np.pi*r*r*np.conj(jn(m,nc*k*r))*jvp(m,nc*k*r))

def timeDelayEigenfunctionIm(r,m,nc,k):
	eps2 = np.imag(nc*nc)
	return np.imag(-1j*k*eps2*2.0*np.pi*r*r*np.conj(jn(m,nc*k*r))*jvp(m,nc*k*r))

# -- We compute the determinant of the S-matrix. 
def sMatrixElementCoefficients(k,n1,n2,m):
  A = -n1*hankel2(m,n2*k)*jvp(m,n1*k)+n2*jn(m,n1*k)*h2vp(m,n2*k)
  B = n1*hankel1(m,n2*k)*jvp(m,n1*k)-n2*jn(m,n1*k)*h1vp(m,n2*k)
  Ap = -n1*n1*jvp(m,n1*k,2)*hankel2(m,n2*k)+n2*n2*jn(m,n1*k)*h2vp(m,n2*k,2)
  Bp = n1*n1*jvp(m,n1*k,2)*hankel1(m,n2*k)-n2*n2*jn(m,n1*k)*h1vp(m,n2*k,2)
  return A, B, Ap, Bp

def sMatrix(k,n1,n2,mMax):
	sMat = np.zeros((2*mMax+1,2*mMax+1))

def qMatrix(k,n1,n2,mMax):
	# -- We construct the scattering matrix and its derivative.
	sMat = np.zeros((2*mMax+1,2*mMax+1), dtype=complex)
	sMati = np.zeros((2*mMax+1,2*mMax+1), dtype=complex)
	sMatp = np.zeros((2*mMax+1,2*mMax+1), dtype=complex)
	eigen = np.zeros((2*mMax+1,2*mMax+1), dtype=complex)
	for i in range(2*mMax+1):
		j = i - mMax
		A, B, Ap, Bp = sMatrixElementCoefficients(k,n1,n2,j)
		sMat[i,i] = A/B
		sMati[i,i] = B/A
		sMatp[i,i] = (Ap*B-A*Bp)/(B*B)
		#eigen[i,i] = quadrature(timeDelayEigenfunctionRe, 0.0, 1.0, (j,n1,k))[0]+1j*1.0*quadrature(timeDelayEigenfunctionIm,0.0,1.0,(j,n1,k))[0]

	return -1j*np.dot(sMati,sMatp)+eigen

# ------------------ MAIN FUNCTION -------------------- #
minK = 9.0
maxK = 9.5
sizeK = 1000
mMax = 25
n1 = 3.2 
n1c = 3.2-1j*0.001	
n2 = 1.0

kValues = np.linspace(minK,maxK,sizeK)
timeDelays = np.zeros((2*mMax+1,sizeK),dtype=complex)
timeDelaysC = np.zeros((2*mMax+1,sizeK),dtype=complex)

for i in range(sizeK):
	# -- Compute the Q-matrix.
	qMat = qMatrix(kValues[i],n1,n2,mMax)
	qMatC = qMatrix(kValues[i]-1j*np.real(n1c)*np.imag(n1c),n1c,n2,mMax)
	timeDelays[:,i] = eigvals(qMat)
	timeDelaysC[:,i] = eigvals(qMatC)

fig = plt.figure(figsize=(8,4))
ax1 = fig.add_subplot(121)
plt.plot(kValues,np.amax(np.real(timeDelays), axis=0))
plt.plot(kValues,np.amax(np.imag(timeDelays), axis=0))
ax1.set_yscale('log')
ax1.set_ylim((1,1e5))

ax2 = fig.add_subplot(122)
plt.plot(kValues,np.amax(np.real(timeDelaysC), axis=0))
plt.plot(kValues,np.amax(np.imag(timeDelaysC), axis=0))
ax2.set_yscale('log')
ax2.set_ylim((1,1e5))
plt.show()