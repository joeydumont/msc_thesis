# -------------------- Information -------------------- #
# Author:       Joey Dumont <joey.dumont@gmail.com>     #
# Date created: Mar. 22nd, 2014							#
# Date mod.:    Mar. 22nd, 2014							#
# Description:  We show the form of the Hankel matrix	#
#				present in the final Hadamard product	#
#				of the RZ algorithm. 					#
# ----------------------------------------------------- #

# --------------- Modules Importation ----------------- #
# Importing numerical analysis package and
# morgenstemning color map.
from pylab import *
from scipy.special import hankel1,hankel2
from matplotlib.colors import LogNorm
import morgenstemning as mrg

# Setting the rc parameters.
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage[charter]{mathdesign}'
rcParams['font.size'] = 10
rcParams['axes.labelsize'] = 'large'
rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['legend.numpoints'] = 3

ms, msi = mrg.morgenstemning()

# --------------- Function Definition ----------------- #
def ratioHankelFunctions(m1, m2, z):
  return hankel2(m2,z)/hankel1(m1,z)


# ---------------- Data Manipulation ------------------ #
# -- Data Structures and Generation
M = 100
hankelMat = zeros((2*M+1,2*M+1), dtype='complex')

for i in range(2*M+1):
  for j in range(2*M+1):
    hankelMat[j,i] = ratioHankelFunctions(-M+j,-M+i,10.0)

# ------------------ Plotting Data -------------------- #
fig1 = figure(figsize=(4,3))
ax1 = fig1.add_subplot(111)

pcolor(abs(hankelMat), cmap=msi, norm=LogNorm(), rasterized=True)
colorbar()

ticks = arange(0, 2*M+1, M/2)
labels = (('$-M$', '$-M/2$', '$0$', '$M/2$', '$M$'))
xticks(ticks,labels)
yticks(ticks,labels[::-1])
xlim((0,2*M))
ylim((0,2*M))
xlabel("$m'$")
ylabel("$m$")
savefig("absHadamard.pdf", bbox_inches='tight')