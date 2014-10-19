# ------------------- Information --------------------- #
# Author:		Joey Dumont <joey.dumont@gmail.com>		#
# Date created:	Jun 28th, 2014							#
# Date mod.		Jun 28th, 2014							#
# Description:	We plot the far-field of the dipole		#
#				antenna. 								#
# ----------------------------------------------------- #

# --------------- Modules Importation ----------------- #
import numpy as np
import matplotlib.pyplot as plt

# Setting the rc parameters.
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = [r"\usepackage[charter]{mathdesign}"]
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 10
plt.rcParams['legend.numpoints'] = 3

# --------------- Function Definition ----------------- #
def weightAlphaTheta(theta, ell):
	return -(np.cos(ell*np.cos(theta))-np.cos(ell))/np.sin(theta)

# ------------------ MAIN FUNCTION -------------------- #
# Half-wave dipole
ell = [np.pi/2]

theta = np.linspace(0.01,2*np.pi,500)

fig = plt.figure(figsize=(3,3))
ax = fig.add_subplot(111, polar=True)
for i in range(len(ell)):
	farField = abs(weightAlphaTheta(theta, ell[i]))**2
	plt.plot(theta,farField)	

ax.set_rmax(1.0)
ax.set_theta_zero_location('N')
plt.annotate("$z$-axis", xy=(0, 0), xytext=(0.18 ,0.8), arrowprops=dict(arrowstyle="<-",facecolor='r'), fontsize=10)

plt.savefig("farField-dipole.pdf", bbox_inches='tight')