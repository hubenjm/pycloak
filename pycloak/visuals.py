import antenna as ant
import controlregion as cr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

def plot_setup(a, c):
	""" (antenna.AntennaArray, controlregion.ControlRegion) -> None
	"""
	for j in range(a.n_antennas):
		x_a = a.get_boundary(j)
		plt.plot(x_a[0,:], x_a[1,:], 'r')
	
	x_cn = c.get_all_near_points() 
	x_cf = c.get_dBR_points()
	plt.plot(x_cn[0,:], x_cn[1,:], 'b')
	plt.plot(x_cf[0,:], x_cf[1,:], 'b')
	plt.tight_layout()
	plt.show()
	

def test1():
	KeyTex = (r'\frac{\| [K \phi]_{1} - f_{1} \|_{L^{2}(\partial D_{c})}}{\|f_{1}\|_{L^{2}(\partial D_{c})}',\
		r"\| [K \phi]_{2} - f_{2} \|_{L^{\infty}(\partial B_{R})} \|",\
		r"\| [K\phi]_{2} - f_{2} \|_{L^{2}(\partial B_{R})}", \
		r"\frac{\| [K \phi]_{1} - f_{1} \|_{L^{2}(\partial D_{c})}}{\|f_{1}\|_{L^{2}(\partial D_{c})} + \| [K \phi]_{2} - f_{2} \|_{L^{\infty}(\partial B_{R})} \|",\
		r"\|\phi\|_{L^{2}(\partial D_{a})}")
	t = np.arange(0., 1.01, 0.01)
	s = np.cos(4*np.pi*t) + 2
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.plot(t,s)
	
	plt.xlabel(r'$$\frac{\| [K \phi]_{1} - f_{1} \|_{L^{2}(\partial D_{c})}}{\|f_{1}\|_{L^{2}(\partial D_{c})}}$$', fontsize=16)
	plt.ylabel(r"$\| [K \phi]_{2} - f_{2} \|_{L^{\infty}(\partial B_{R})}$", fontsize=16)
	plt.title(r"$\|\phi\|_{L^{2}(\partial D_{a})}$", fontsize=16)
	
	plt.tight_layout()
	
	plt.savefig('test.png')
	plt.show()
	
def test2():
	a = ant.AntennaArray(ant.circle_s, ant.circle_ds, 0.3, 5, 64)
	c = cr.ControlRegion(0.11, 0.2, -np.pi/4, np.pi/4, 0, 10, 512,128,128,16)
	plot_setup(a,c)
	

if __name__=="__main__":
	print "Running tests..."	
	test2()	
