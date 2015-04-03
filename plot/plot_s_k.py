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
	
	NFE = near field relative error
	FFE = far field average L2 error
	F = value of mismatch functional
	PHINORM = norm of phi
	FLAG = indicates how Newton's method terminated
	ITER = number of Newton's method iterations
	ALPHA = value of optimal regularization parameter
	
	"""
	data = np.load(filename)
	Alpha = data['ALPHA']
	Cn = data['NFE']
	Cf = data['FFE']
	P = data['PHINORM']
	Flag = data['FLAG']
	Iter = data['ITER']
	F = data['F']
	K = data['K']
	S = data['S']
	delta = data['DELTA']

	#auxiliary variables
	LogA = np.log10(Alpha)
	LogP = np.log10(P)

	plt.rc('text', usetex=False)
	plt.rc('font', family='serif')
	orient = 'vertical'
	XLABEL = r"$s$"
	YLABEL = r"$k$"
	
	Ns = S.shape[0]
	Nk = S.shape[1]
	N = Ns*Nk
	
	S = np.reshape(S, (N, ))
	K = np.reshape(K, (N, ))
	LogA = np.reshape(LogA, (N, ))
	LogP = np.reshape(LogP, (N, ))
	Cn = np.reshape(Cn, (N, ))
	Cf = np.reshape(Cf, (N, ))
	

	fig = plt.figure(1)
	ax = fig.add_subplot(221, projection='3d')
	ax.set_xlabel(XLABEL, fontsize=AXESFONTSIZE)
	ax.set_ylabel(YLABEL, fontsize=AXESFONTSIZE)
	ax.set_zlabel(r"$\log_{10}{\alpha}$", fontsize=AXESFONTSIZE)
	ax.set_xlim3d(S.min(), S.max())
	ax.set_ylim3d(K.min(), K.max())
	ax.set_zlim3d(LogA.min(), LogA.max())
	surf1 = ax.plot_trisurf(S, K, LogA, vmin=LogA.min(), vmax=LogA.max(), cmap=cm.jet, edgecolor='none')
	#surf1 = ax.plot_surface(S, K, LogA, vmin=LogA.min(), vmax=LogA.max(), cmap=cm.jet, edgecolor='none')
	fig.colorbar(surf1, shrink=0.5, aspect=5, pad=0.05, orientation = orient)

	ax = fig.add_subplot(222, projection='3d')
	ax.set_xlabel(XLABEL, fontsize=AXESFONTSIZE)
	ax.set_ylabel(YLABEL, fontsize=AXESFONTSIZE)
	ax.set_zlabel("\n" + r"$\log\left(\Vert \phi\Vert_{L^{2}(\partial D_{a})} \right)$", fontsize=AXESFONTSIZE, linespacing=2)
	ax.set_xlim3d(S.min(), S.max())
	ax.set_ylim3d(K.min(), K.max())
	surf2 = ax.plot_trisurf(S, K, LogP, cmap=cm.jet, vmin=LogP.min(), vmax=LogP.max(), edgecolor='none')
	#surf2 = ax.plot_surface(S, K, LogP, cmap=cm.jet, vmin=LogP.min(), vmax=LogP.max(), edgecolor='none')
	fig.colorbar(surf2, shrink=0.5, aspect=5, pad=0.05, orientation = orient)

	ax = fig.add_subplot(223, projection='3d')
	ax.set_xlabel(XLABEL, fontsize=AXESFONTSIZE)
	ax.set_ylabel(YLABEL, fontsize=AXESFONTSIZE)
	ax.set_zlabel("\n" + r"$\frac{\Vert K_{1}\phi - f_{1}\Vert}{\Vert f_{1}\Vert}$", fontsize=AXESFONTSIZE, linespacing=2)
	ax.set_xlim3d(S.min(), S.max())
	ax.set_ylim3d(K.min(), K.max())
	surf3 = ax.scatter(S, K, Cn, c=Cn, s = 5)
	#surf3 = ax.plot_surface(S, K, Cn, cmap=cm.jet, vmin=Cn.min(), vmax=Cn.max(), edgecolor='none')
	fig.colorbar(surf3, shrink=0.5, aspect=5, pad=0.05, orientation = orient)

	ax = fig.add_subplot(224, projection='3d')
	ax.set_xlabel(XLABEL, fontsize=AXESFONTSIZE)
	ax.set_ylabel(YLABEL, fontsize=AXESFONTSIZE)
	ax.set_zlabel("\n" + r"$\frac{1}{\sqrt{2\pi R}}\log{\Vert K_{2}\phi \Vert}$", fontsize=AXESFONTSIZE, linespacing=2)
	ax.set_xlim3d(S.min(), S.max())
	ax.set_ylim3d(K.min(), K.max())
	surf4 = ax.plot_trisurf(S, K, Cf, cmap=cm.jet, vmin=Cf.min(), vmax=Cf.max(), edgecolor='none')
	#surf4 = ax.plot_surface(S, K, Cf, cmap=cm.jet, vmin=Cf.min(), vmax=Cf.max(), edgecolor='none')
	fig.colorbar(surf4, shrink=0.5, aspect=5, pad=0.05, orientation = orient)
	ax.view_init(azim=-129, elev=20)

	fig.subplots_adjust(left = 0.02, right = 0.98, wspace=0.05, hspace=0.05, top=0.95, bottom=0.05)
	plt.show()
