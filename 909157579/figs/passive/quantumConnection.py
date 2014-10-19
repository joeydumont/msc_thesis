# ------------------- Information --------------------- #
# Author:		Joey Dumont <joey.dumont@gmail.com>		#
# Date created:	Jun. 23rd, 2014							#
# Date mod.	:	Jun. 23rd, 2014							#
# Description:	We plot the potential of the Helmholtz-	#
#				Schrödinger comparison. 				#
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
def potential(r,k,m,nc):
	if r > 1:
		pot = (m*m-0.25)/(r*r)
	else:
		pot = k*k*(1-nc*nc)+(m*m-0.25)/(r*r)

	return pot

fig = plt.figure(figsize=(3,3))
ax = fig.add_subplot(111)

size = 1000
r1 = np.linspace(0.01,1,size)
r2 = np.linspace(1.01,3,size)
pot1 = np.zeros((size))
pot2 = np.zeros((size))
for i in range(size):
	pot1[i] = potential(r1[i], 5.0, 10, 2)
	pot2[i] = potential(r2[i], 5.0, 10, 2)

plotArgs = {"linestyle": '-', 'color': 'k', 'linewidth': 1.5}
plt.plot(r1[:],pot1[:], **plotArgs)
plt.fill_between(r1[:],pot1[:],0, color='0.85')
plt.plot(r2[:],pot2[:], **plotArgs)
plt.fill_between(r2[:],pot2[:],0, color='0.85')

ax.set_xlabel(r"$r/R_0$")
ax.set_ylim((0,150))
ax.yaxis.set_visible(False)
plt.axvline(x=1, ymin=25/150, ymax=100/155, color='k', linestyle='--', linewidth=1)

plt.savefig("quantumPotential.pdf", bbox_inches='tight')