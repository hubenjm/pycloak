import matplotlib.pyplot as plt
from matplotlib import cm, ticker
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import sys

AXESFONTSIZE = 24
TEXTFONTSIZE = 26
TITLEFONTSIZE = 34
TICKFONTSIZE = 22
LineWidth = 0.3
    
def plot(filename):
	"""
	The data is stored in a numpy.npz file which is easily reloadable directly to numpy arrays.
	
	D = distance from control region to antenna
	EPSILON = relative noise level added to data on D_c
	NFE = near field relative error
	FFE = far field average L2 error
	F = value of mismatch functional
	PHINORM = norm of phi
	STABILITY = ||phi_epsilon - phi_0||/||phi_0||
	FLAG = indicates how Newton's method terminated
	ITER = number of Newton's method iterations
	ALPHA = value of optimal regularization parameter
	
	"""
	data = np.load(filename)
	Alpha = data['ALPHA']
	Cn = data['NFE']
	Cf = data['FFE']
	P = data['PHINORM']
	Sr = data['STABILITY']
	Flag = data['FLAG']
	Iter = data['ITER']
	F = data['F']
	D = data['D']
	E = data['EPSILON']
	delta = data['DELTA']

	#auxiliary variables
	LogA = np.log10(Alpha)
	LogP = np.log10(P)

	plt.rc('text', usetex=False)
	plt.rc('font', family='serif')
	orient = 'vertical'
	XLABEL = r"$d$"
	YLABEL = r"$\epsilon$"
	
	Nd = D.shape[0]
	Nepsilon = D.shape[1]
	N = Nd*Nepsilon
	
	D = np.reshape(D, (N, ))
	E = np.reshape(E, (N, ))
	Sr = np.reshape(Sr, (N, ))
	LogA = np.reshape(LogA, (N, ))
	LogP = np.reshape(LogP, (N, ))
	G = np.sqrt(np.abs(F))/delta
	G = np.reshape(G, (N, ))

	fig = plt.figure(1)
	ax = fig.add_subplot(221, projection='3d')
	ax.set_xlabel(XLABEL, fontsize=AXESFONTSIZE)
	ax.set_ylabel(YLABEL, fontsize=AXESFONTSIZE)
	ax.set_zlabel(r"$\log_{10}{\alpha}$", fontsize=AXESFONTSIZE)
	ax.set_xlim3d(D.min(), D.max())
	ax.set_ylim3d(E.min(), E.max())
	ax.set_zlim3d(LogA.min(), LogA.max())
	surf1 = ax.plot_trisurf(D, E, LogA, vmin=LogA.min(), vmax=LogA.max(), cmap=cm.jet, edgecolor='none')
	fig.colorbar(surf1, shrink=0.5, aspect=5, pad=0.05, orientation = orient)

	ax = fig.add_subplot(222, projection='3d')
	ax.set_xlabel(XLABEL, fontsize=AXESFONTSIZE)
	ax.set_ylabel(YLABEL, fontsize=AXESFONTSIZE)
	ax.set_zlabel("\n" + r"$\log\left(\Vert \phi\Vert_{L^{2}(\partial D_{a})} \right)$", fontsize=AXESFONTSIZE, linespacing=2)
	ax.set_xlim3d(D.min(), D.max())
	ax.set_ylim3d(E.min(), E.max())
	surf2 = ax.plot_trisurf(D, E, LogP, cmap=cm.jet, vmin=LogP.min(), vmax=LogP.max(), edgecolor='none')
	fig.colorbar(surf2, shrink=0.5, aspect=5, pad=0.05, orientation = orient)

	ax = fig.add_subplot(223, projection='3d')
	ax.set_xlabel(XLABEL, fontsize=AXESFONTSIZE)
	ax.set_ylabel(YLABEL, fontsize=AXESFONTSIZE)
	ax.set_zlabel("\n" + r"$\frac{\Vert \phi^{\epsilon} - \phi_{0}\Vert}{\Vert \phi^{0}\Vert}$", fontsize=AXESFONTSIZE, linespacing=2)
	ax.set_xlim3d(D.min(), D.max())
	ax.set_ylim3d(E.min(), E.max())
	surf3 = ax.plot_trisurf(D, E, Sr, cmap=cm.jet, vmin=Sr.min(), vmax=Sr.max(), edgecolor='none')
	fig.colorbar(surf3, shrink=0.5, aspect=5, pad=0.05, orientation = orient)

	ax = fig.add_subplot(224, projection='3d')
	ax.set_xlabel(XLABEL, fontsize=AXESFONTSIZE)
	ax.set_ylabel(YLABEL, fontsize=AXESFONTSIZE)
	ax.set_zlabel("\n" + r"$\frac{\sqrt{|F(\alpha)|}}{\delta}$", fontsize=AXESFONTSIZE, linespacing=2)
	ax.set_xlim3d(D.min(), D.max())
	ax.set_ylim3d(E.min(), E.max())
	scatter4 = ax.scatter(D, E, G, c=G, s=5)
	fig.colorbar(scatter4, shrink=0.5, aspect=5, pad=0.05, orientation = orient)
	ax.view_init(azim=-129, elev=20)

	fig.subplots_adjust(left = 0.02, right = 0.98, wspace=0.05, hspace=0.05, top=0.95, bottom=0.05)
	plt.show()
