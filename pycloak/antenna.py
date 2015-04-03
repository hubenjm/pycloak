from math import pi
import numpy as np

class Antenna:
	"""
	General super class for all antennas to inherit from
	"""
	def integrate(self, A, dim=0):
		"""
		Integrate the array A over the antenna along dimension dim
		"""
		assert False, 'method must be defined!'
		
	def norm(self, f, dim=0):
		""" (numpy.array) -> float
		Takes a vector phi of appropriate length and computes the L2 norm of phi as a function on the antenna boundary.
		"""
		return np.sqrt(self.integrate(f*np.conjugate(f)))
		
	def nu(self):
		"""
		return a 2 by self.n_points() array of the unit normal vectors to boundary of Antenna
		"""
		assert False
		
	def get_boundary(self):
		"""
		return a 2 by self.n_points() array of points along the boundary of the Antenna
		"""
		assert False
		
	def dl(self):
		"""
		return a self.n_points() length array of arc length differentials
		"""
		assert False
		
	def n_points(self):
		assert False


class Rectangle(Antenna):
	"""Consists of a single rectangular shaped antenna, centered at (0,0)
	"""
	def __init__(self, width, height, n, xc = np.array([0,0])):
		self.width = width
		self.height = height
		self.xc = xc #center coordinates of antenna
		self.n = n #number of sample points around whole rectangle
		
	def perimeter(self):
		return 2.*self.width + 2.*self.height
	
	def ds(self):
		return self.perimeter()/self.n_points()
	
	def dl(self):
		return self.ds()*np.ones(self.n_points())	

	def get_boundary(self):
		"""Returns a 2 by self.n_points() matrix whose columns are evenly spaced points along the antenna with respect to arclength.
		"""
		L = np.linspace(0, self.perimeter() - self.ds(), self.n)
		j1 = len(L[L<self.height])
		j2 = len(L[L<self.height + self.width])
		j3 = len(L[L<2*self.height + self.width])
		j4 = self.n

		a1 = -self.width/2.0*np.ones(j1)
		a2 = self.width/2.0*np.ones(j3-j2)
		b1 = self.height/2.0*np.ones(j2-j1)
		b2 = -self.height/2.0*np.ones(j4-j3)

		X = np.hstack([a1, L[j1:j2]-self.height - self.width/2.0, a2, -L[j3:]+2*self.height+1.5*self.width])
		Y = np.hstack([L[:j1]-self.height/2.0, b1, -L[j2:j3]+1.5*self.height + self.width, b2])
		return np.vstack([X,Y]) + np.reshape(self.xc, (2,1))
		
	def integrate(self, A, dim=0):
		"""
		Integrate the matrix A on the antenna along dimension dim
		"""
		m = len(A.shape)
		assert m <= 2 and m > 0
		assert dim in range(m) and A.shape[dim] == self.n_points()
		
		return np.sum(self.ds()*A, dim)
	
	def nu(self):
		L = np.linspace(0, self.perimeter() - self.ds(), self.n)
		j1 = len(L[L<self.height])
		j2 = len(L[L<self.height + self.width])
		j3 = len(L[L<2*self.height + self.width])
		j4 = self.n
		
		nu1 = np.vstack([-1,0])
		nu2 = np.vstack([0,1])
		nu3 = np.vstack([1,0])
		nu4 = np.vstack([0,-1])
		
		return np.hstack([np.tile(nu1, j1), np.tile(nu2, j2-j1), np.tile(nu3, j3-j2), np.tile(nu4, j4-j3)])
		
	def n_points(self):
		return self.n
	
class PolarArray(Antenna):
	"""Consists of a single antenna profile given by s repeated 2n+1 times with spacing <spacing> and with the middle antenna centered at the origin. We index the antennas in the array starting at j=-n and going up to j=n. Thus j=0 corresponds to the antenna centered at (0,0)
	"""
	ptolerance = 1e-12 #tolerance for how close to periodic s(theta) needs to be
	def __init__(self, s, ds, centers, n):
		assert abs(s(2*pi) - s(0)) < self.ptolerance and abs(ds(2*pi) - ds(0)) < self.ptolerance, "Inputs s and ds are not 2*pi-periodic"
		self.s = s
		self.ds = ds #note that s and ds must be defined for np.arrays
		self.n = n #parameter to determine how many angle samples on each antenna boundary
					#also determines how many basis functions to use on each antenna
		if centers.shape == (2,):
			#only one antenna in array
			self.spacing = 0
			self.n_antennas = 1
			self.centers = np.reshape(centers, (2,1))
		else:
			self.spacing = np.sqrt(np.sum((centers[:,1:] - centers[:,:-1])**2, 0)) #vector of spacings
			self.n_antennas = centers.shape[1]
			self.centers = centers
	
	def n_points(self):
		return self.n*self.n_antennas
	
	def get_center(self, j):
		"""int -> numpy.array
		Takes an int j between 0 and self.n_antennas-1 as input and returns the shape (2,1) array
		of the coordinates of the jth antenna center.
		"""
		return self.centers[:,j]
			
	def get_centers(self):
		"""returns a 2 by self.n_antennas array of the center of each antenna in the array
		"""
		return self.centers
	
	def get_thetas(self):
		""" -> numpy.array
		Returns an array of n regularly spaced angles on the interval [0,2*pi]
		"""
		return np.linspace(0, 2*np.pi*(1. - 1./self.n), self.n)
	
	def dtheta(self):
		return 2*np.pi/self.n
		
	def eval_s(self):
		"""Here thetas is a numpy array of angles (or a float) 
		"""
		return self.s(self.get_thetas())
		
	def eval_ds(self):
		"""Evaluate the derivative of s (given by self.ds) at the points in the numpy array thetas (thetas can also be a single number)
		"""
		return self.ds(self.get_thetas())
		
	def dl(self, part=False):
		"""
		Differential on boundary. Should be a self.n_points() by 1 np.array. If part=True, then returns
		an array of length self.n (i.e. when dealing with just one antenna in array)
		"""
		if part:
			return self.dtheta()*np.sqrt(self.eval_ds()**2 + self.eval_s()**2)
		else:
			return self.dtheta()*np.tile(np.sqrt(self.eval_ds()**2 + self.eval_s()**2), self.n_antennas)
		
	def nu(self):
		return np.tile(self.nu_j(), self.n_antennas)
	
	def nu_j(self):
		thetas = self.get_thetas()
		S = self.eval_s()
		Sprime = self.eval_ds()
		J = np.sqrt(S**2 + Sprime**2) #shape (self.n,)
		return np.array([S*np.cos(thetas) + Sprime*np.sin(thetas), S*np.sin(thetas) - Sprime*np.cos(thetas)])/np.array([J, J]) #nu.shape = (2,n)	
		
	def get_boundary_j(self, j):
		"""Given thetas is a vector of length n, then we return a 2 by n matrix whose columns are the coordinates of the boundary of the jth antenna in the array evaluated at the angles given in thetas.
		"""
		assert j >= 0 and j < self.n_antennas
		S = self.eval_s()
		thetas = self.get_thetas()
		return(np.array([S*np.cos(thetas), S*np.sin(thetas)]) + np.reshape(self.get_center(j), (2,1)))
		 
	def get_boundary(self):
		return np.hstack([self.get_boundary_j(j) for j in range(self.n_antennas)])
		
	def integrate(self, A, dim=0):
		"""
		Integrate the matrix A on the antenna along dimension dim. A should have shape (D_a.n_points(),) or (N, D_a.n_points())
		"""
		m = len(A.shape)
		assert m <= 2 and m > 0
		assert dim in range(m) and A.shape[dim] == self.n_points()

		if dim == 0:
			return np.dot(self.dl(), A)
		else:
			return np.dot(A, self.dl())
	
class Polar(PolarArray):
	"""This is a class to represent a given Polar object by a boundary defining function f(theta) and its derivative f'(theta). We are restricting ourselves here to antenna which have a polar function representation. We assume (but don't check) that s and ds are periodic on [0, 2*pi], with numpy arrays as inputs"""
	def __init__(self, s, ds, center, n):
		#center has shape (2,)
		PolarArray.__init__(self, s, ds, center, n)

def make_circle_s(radius):
	def f(theta): 
		theta = np.asarray(theta)
		return radius*np.ones(theta.shape)
	return f

def circle_ds(theta):
	theta = np.asarray(theta)
	return np.zeros(theta.shape)
	
def make_centers(x1, x2, t0, t1, n):
	"""
	Here, x1 and x2 are a function of t
	"""
	t = np.linspace(t0,t1,n)
	return np.vstack([x1(t), x2(t)])

def test1():
	radius = 0.1
	A = Polar(make_circle_s(radius), circle_ds, np.array([0,0]), 128)
	thetas = A.get_thetas()
	assert len(thetas) == A.n
	assert A.n_antennas == 1
	flist = [make_circle_s(radius), circle_ds]
	for f in flist:
		print f(thetas)
	print A.eval_s()
	print A.nu()
	
	print A.ptolerance
	B = Polar(lambda x: 0.1*np.ones(np.asarray(x).shape), lambda x: np.zeros(np.asarray(x).shape), np.array([0,0]), 16)
	print B.eval_s()
	print B.eval_ds()
	
	C = PolarArray(make_circle_s(radius), circle_ds, make_centers(np.cos, np.sin, 0, np.pi, 10), 128)
	print C.get_boundary_j(3).shape
	print C.get_centers().shape
	print isinstance(C, Polar)
	print isinstance(A, Polar)

def test2():
	import matplotlib.pyplot as plt
	n = 128
	w = 1
	h = 3
	print("Testing rectangular antenna...")
	D = Rectangle(w,h,n, np.array([1,2]))
	#print D.get_boundary()
	phi = np.random.rand(n)
	print(D.norm(phi))
	assert D.norm(phi) == np.sqrt(np.sum(phi*np.conjugate(phi)*D.ds()))
	fig = plt.figure(1)
	ax = fig.add_subplot(111)
	points = D.get_boundary()
	print D.nu()
	ax.scatter(points[0,:], points[1,:], s=0.3)
	ax.set_xlim(-10,10)
	ax.set_ylim(-10,10)
	ax.grid()
	plt.show()

if __name__=="__main__":
	print "Running tests..."	
	test1()
	test2()
