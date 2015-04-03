import numpy as np
import controlregion as cr
from scipy.special import hankel1

class Source:
	def trace(self, D_c):
		assert False

class PointSource(Source):
	def __init__(self, x0):
		assert x0.size == 2
		self.x0 = x0
		
	def trace(self, D_c, k):
		""" (controlregion.ControlRegion) -> numpy.array
		"""
		x = D_c.get_boundary()
		return (1j/4)*hankel1(0,k*np.linalg.norm(x - np.reshape(self.x0, (2,1)), axis=0))

class PlaneWave(Source):
	def __init__(self, xi):
		self.xi = xi
		self.xihat = np.array([np.cos(xi), np.sin(xi)])
	
	def trace(self, D_c, k):
		x = D_c.get_boundary()
		return np.exp(1j*k*(self.xihat[0]*x[0,:] + self.xihat[1]*x[1,:]))

def test1():
	k = 2
	x0 = np.array([10000,0])
	x = np.random.rand(2,1000)	
	xi = np.pi/2
	p = PlaneWave(xi)
	param = cr.make_annularsector_param_dict(1.11,  1.14, -np.pi/4, np.pi/4, 0, 256, 256, 64)
	D_c = cr.AnnularSector(param)
	vals = p.trace(D_c, k)
	print vals.shape
	print vals

if __name__=="__main__":
	print "Running tests..."	
	test1()		
