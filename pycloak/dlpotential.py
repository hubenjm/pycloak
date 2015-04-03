"""Set of routines to compute 
"""
import antenna
import controlregion
import numpy as np
import boundaryfunction
from scipy.special import hankel1
from scipy.integrate import simps

class ControlProblem:
	"""
	"""
	def __init__(self, D_a, D_c, D_far, u0, k, basis_type='delta'):
		assert k >= 0, 'k must be non-negative'
		assert basis_type in ('delta', 'fourier'), 'basis_type must be either fourier or delta'
		if basis_type == 'fourier':
			assert isinstance(D_a, (antenna.PolarArray, antenna.Polar))
		assert isinstance(u0, boundaryfunction.Source)	
			
		self.D_a = D_a
		self.D_c = D_c
		self.D_far = D_far
		self.u0 = u0 #function to match on D_c
		self.basis_type = basis_type
		self.k = k
		self.A_near = initialize_A_near(k, D_a, D_c, basis_type)
		self.A_far = initialize_A_far(k, D_a, D_far, basis_type)
		
		self.A_near_iscurrent = True
		self.A_far_iscurrent = True
		
		#store A_near and A_far in memory to avoid recomputing 
	
	def set_k(self, k):
		self.k = k
		self.A_near_iscurrent = False
		self.A_far_iscurrent = False
		
	def set_D_c_parameters(self, param):
		"""
		We assume param is a dictionary object
		"""
		for j in param:
			if j in self.D_c.param:
				self.D_c.param['j'] = param['j']
		self.A_near_iscurrent = False
		
	def set_D_far_parameters(self, param):
		for j in param:
			if j in self.D_far.param:
				self.D_fac.param['j'] = param['j']
		self.A_far_iscurrent = False
	
	def A(self):
		return np.concatenate((self.compute_A_near(), self.compute_A_far()), axis=0)
			
	def B_near(self):
		return np.dot(np.conjugate(self.compute_A_near()).T, self.compute_A_near())
	
	def B_far(self):
		return np.dot(np.conjugate(self.compute_A_far()).T, self.compute_A_far())
		
	def B(self):
		return self.B_near() + self.B_far()
	
	def error(self, g, f):
		N = self.D_c.n_points() + self.D_far.n_points()
		assert len(f) == N and len(g) == N
		nc = self.D_c.n_points()
		residual = g - f
		near_error_relative = self.D_c.L2(residual[:nc])/self.D_c.L2(f[:nc])
		far_error_l2ave = self.D_far.L2_ave(residual[nc:])
		
		return near_error_relative + far_error_l2ave
	
	def compute_A_near(self):
		""" 
		Computes the near matrix representation of the DL operator K with respect
		to either the product space of exp basis functions exp(i*l*theta) for l=0..n-1,
		or a delta/spline basis, depending on the value of self.basis_type. The return object
		is a numpy.array of size self.control_space.D_c.n_points() by self.D_a.n_points().
	
		C.k = wave number. If C.k > 0 then we use the Helmholtz fundamental solution.
		Otherwise, if C.k = 0 we use the fundamental solution for the Laplace equation.
		"""
		if self.A_near_iscurrent:
			return self.A_near
			print "A_near is current"
		else:
			nc = self.D_c.n_points()
			na = self.D_a.n_points()
			A_near = np.zeros((nc, na), dtype=complex)
			X = self.D_c.get_boundary()
	
			if self.basis_type == 'delta':
				for j in range(nc):
					x_j = X[:,j]
					A_near[j,:] = self.D_a.dl()*ddn_helmholtz_fs(self.k, self.D_a, x_j)
		
			if self.basis_type == 'fourier':
				for i in range(self.D_a.n_antennas):
					for j in range(nc):
						x_j = X[:,j]
						A_near[j,i*self.D_a.n:(i+1)*self.D_a.n] = self.D_a.n*np.fft.ifft(ddn_helmholtz_fs_j(self.k, self.D_a, x_j, i)*self.D_a.dl(part=True))	
			self.A_near = A_near
			self.A_near_iscurrent = True	
			return self.A_near
			
			
	def set_D_a(self, D_a_new):
		self.D_a = D_a_new
		self.A_near_iscurrent = False
		self.A_far_iscurrent = False
	
	def compute_A_far(self):
		if self.A_far_iscurrent:
			return self.A_far
		else:
			nfar = self.D_far.n_points()
			na = self.D_a.n_points()
			A_far = np.zeros((nfar, na), dtype=complex)
			X = self.D_far.get_boundary()
	
			if self.basis_type == 'delta':
				for j in range(nfar):
					x_j = X[:,j]
					A_far[j,:] = self.D_a.dl()*ddn_helmholtz_fs(self.k, self.D_a, x_j)
		
			if self.basis_type == 'fourier':
				for i in range(self.D_a.n_antennas):
					for j in range(nfar):
						x_j = X[:,j]
						A_far[j,i*self.D_a.n:(i+1)*self.D_a.n] = self.D_a.n*np.fft.ifft(ddn_helmholtz_fs_j(self.k, self.D_a, x_j, i)*self.D_a.dl(part=True))	
			self.A_far = A_far
			self.A_far_iscurrent = True
			return self.A_far
		
	def near_trace(self, phi):
		""" numpy.array -> numpy.array
		Here we assume phi is an np.array of shape (self.D_a.n_points(), ) that represents
		its function values on the boundary of the antenna
		"""
		assert len(phi.shape) == 1 and len(phi) == self.D_a.n_points()
	
		na = self.D_a.n_points()
		nc = self.D_c.n_points()
		f_near = np.zeros(nc, dtype=complex)
		X = self.D_c.get_boundary()
		
		if self.k > 0:
			for j in range(nc):
				x_j = X[:,j]
				f_near[j] = self.D_a.integrate(phi*ddn_helmholtz_fs(self.k, self.D_a, x_j), 0)
		if self.k == 0:
			return False #add Laplace later
		
		return f_near
		
	def far_trace(self, phi):
		""" numpy.array -> numpy.array
		"""
		assert len(phi.shape) == 1 and len(phi) == self.D_a.n_points()
		na = len(phi)
		nfar = self.D_far.n_points()
		f_far = np.zeros(nfar, dtype=complex)
		X = self.D_far.get_boundary()
	
		if self.k > 0:
			for j in range(nfar):
				x_j = X[:,j]
				f_far[j] = self.D_a.integrate(phi*ddn_helmholtz_fs(self.k, self.D_a, x_j), 0)
			
		if self.k == 0:
			return False #add Laplace later
		
		return f_far
		
	def control_trace(self, phi):
		f_near = self.near_trace(phi)
		f_far = self.far_trace(phi)
		return np.concatenate((f_near, f_far), axis=0)
		
	def nfe(self, h, f_near):
		nc = self.D_c.n_points()
		na = self.D_a.n_points()
		assert len(f_near) == nc
		assert len(h) == na
		
		A_near = self.compute_A_near()
		residual = np.dot(A_near, h) - f_near
		return np.sqrt(self.D_c.L2(residual)/self.D_c.L2(f_near))
		
	def ffe(self, h, f_far):
		nfar = self.D_far.n_points()
		na = self.D_a.n_points()
		assert len(f_far) == nfar
		assert len(h) == na
		
		A_far = self.compute_A_far()
		residual = np.dot(A_far, h) - f_far
		return np.sqrt(self.D_far.L2_ave(residual))
		
	def F(self, h, f, delta):
		nc = self.D_c.n_points()
		nfar = self.D_far.n_points()
		na = self.D_a.n_points()
		assert len(f) == nc + nfar
		assert len(h) == na
		
		f_near = f[:nc]
		f_near_normsq = self.D_c.L2(f_near)
		
		residual = np.dot(self.A(), h) - f
		s = self.D_c.L2(residual[:nc])/f_near_normsq + self.D_far.L2_ave(residual[nc:]) - delta**2
		return np.real(s)
			
	def invert_coeff_vector(self, w):
		""" numpy.array -> numpy.array
	
		w is the vector of coefficients of a function defined on self.D_a with respect to the basis specified in self. If self.basis_type == 'fourier', we apply the ifft to each part of w to return the actual antenna profile. If self.basis_type == 'delta', we don't need to do anything; just return w as is.
		"""
		assert len(w.shape) == 1 and len(w) == self.D_a.n_points()
		na = len(w)
		if self.basis_type == 'delta':
			return w
	
		if self.basis_type == 'fourier':
			w = np.reshape(w, (self.D_a.n_antennas, self.D_a.n))
			return self.D_a.n*np.reshape(np.fft.ifft(w), na)
	
	def compute_coeff_vector(self, phi):
		assert len(phi.shape) == 1 and len(phi) == D_a.n_points()
		na = len(phi)
	
		if basis_type == 'delta':
			return phi
	
		if basis_type == 'fourier':
			phi = np.reshape(phi, (self.D_a.n_antennas, self.D_a.n))
			return (1./self.D_a.n)*np.reshape(np.fft.fft(phi, axis=1), na)
		
#########################################################################		

def initialize_A_near(k, D_a, D_c, basis_type):
	nc = D_c.n_points()
	na = D_a.n_points()
	A_near = np.zeros((nc, na), dtype=complex)
	X = D_c.get_boundary()

	if basis_type == 'delta':
		for j in range(nc):
			x_j = X[:,j]
			A_near[j,:] = D_a.dl()*ddn_helmholtz_fs(k, D_a, x_j)
	
	if basis_type == 'fourier':
		assert isinstance(D_a, (antenna.Polar, antenna.PolarArray))
		for i in range(D_a.n_antennas):
			for j in range(nc):
				x_j = X[:,j]
				A_near[j,i*D_a.n:(i+1)*D_a.n] = D_a.n*np.fft.ifft(ddn_helmholtz_fs_j(k, D_a, x_j, i)*D_a.dl(part=True))
				
	return A_near
		
def initialize_A_far(k, D_a, D_far, basis_type):
	nfar = D_far.n_points()
	na = D_a.n_points()
	A_far = np.zeros((nfar, na), dtype=complex)
	X = D_far.get_boundary()

	if basis_type == 'delta':
		for j in range(nfar):
			x_j = X[:,j]
			A_far[j,:] = D_a.dl()*ddn_helmholtz_fs(k, D_a, x_j)
	
	if basis_type == 'fourier':
		for i in range(D_a.n_antennas):
			for j in range(nfar):
				x_j = X[:,j]
				A_far[j,i*D_a.n:(i+1)*D_a.n] = D_a.n*np.fft.ifft(ddn_helmholtz_fs_j(k, D_a, x_j, i)*D_a.dl(part=True))
				
	return A_far		


def helmholtz_fs(k, x, y):
	""" (float, numpy.array, numpy.array) -> numpy.array
	
	Compute the value of the fundamental solution to the Helmholtz equation with wave number k. 
	x is a 2 by D_c.n_points() array of points along a particular ControlRegion object. 
	y is a 2 by D_a.n_points() array of points along the boundary of an Antenna object. 
	Output is an D_a.n_points() by D_c.n_points() array of values. 
	The columns correspond to a fixed point in the control region (x) while varying the point on the antenna array (y)
	"""
	assert x.shape[0]==2 and y.shape[0]==2 and k > 0
	nc = x.shape[1]
	na = y.shape[1]
	result = np.zeros((nc,na), dtype=complex)
	
	for i in range(nc):
		result[i,:] = 1j/4*hankel1(0,k*np.linalg.norm(np.reshape(x[:,i],(2,1))-y, axis=0))
	
	return result.T #shape is (na, nc)
	
def ddn_helmholtz_fs(k, D_a, x):
	""" (float, Antenna, numpy.array) -> numpy.array
	
	We assume that x is a single 2-vector not on the Antenna D_a.
	Returns a D_a.n_points() by 1 array (shape = (D_a.n_points(), ))representing the values dPhi/dn(y,x) for each value y on the antenna array D_a.
	"""
	assert x.size==2 and k > 0
	x = np.reshape(x,2) #makes sure that x has shape (2,)
	
	#need to compute normal vectors to each boundary of antenna array
	y = D_a.get_boundary() #y.shape = (2, D_a.n_points())
	nu = D_a.nu() #shape = (2, D_a.n_points())
	
	ynu = np.sum(y*nu, 0) #Ynu.shape = (D_a.n_points(),)
	xnu = np.dot(nu.T, x) #xnu.shape = (D_a.n_points(),)
	nu_diff_mat = ynu - xnu #nu_diff_mat.shape = (D_a.n_points(),)
	dist_array = np.linalg.norm(y-np.reshape(x,(2,1)), axis=0) #dist_array.shape = (D_a.n_points(),)
	return -(1j*k/4)*nu_diff_mat/dist_array*hankel1(1,k*dist_array)

def ddn_helmholtz_fs_j(k, D_a, x, j):
	""" (number, PolarArray, numpy.array, int) -> numpy.array
	
	This is a more compartmentalized version of ddn_helmholtz_fs (works only for PolarAntennaArray).
	Here we just want to output the values of the normal derivative for y on the boundary of
	the jth antenna where j ranges from 0 to a.n_antennas-1.
	
	Again we expect that x is a 2-vector (i.e. x.size = 2). k is also positive.
	"""
	assert isinstance(D_a, (antenna.PolarArray, antenna.Polar))
	assert x.size==2 and k > 0
	x = np.reshape(x,2) #makes sure that x has shape (2,)
	nu = D_a.nu_j()
	y = D_a.get_boundary_j(j)
		
	ynu = np.sum(y*nu, 0) #Ynu.shape = (n,)
	xnu = np.dot(nu.T, x) #xnu.shape = (n,)
	nu_diff_mat = ynu - xnu #nu_diff_mat.shape = (n,)
	dist_array = np.linalg.norm(y-np.reshape(x,(2,1)), axis=0) #dist_array.shape = (n,)
	
	return -(1j*k/4)*nu_diff_mat/dist_array*hankel1(1,k*dist_array)
	
def test1():
	x1 = np.array([2,1])
	x2 = np.array([[2,1]])
	k = 1
	radius = 0.1
	D_a = antenna.PolarArray(antenna.make_circle_s(radius), antenna.circle_ds, antenna.make_centers(np.cos, np.sin, 0, np.pi/2, 5), 32)
	
	centers = D_a.get_centers()
	print "centers.shape =", centers.shape
	print "centers[:,0] =", centers[:,0]
	
	w1 = ddn_helmholtz_fs(k, D_a, x2)
	print "w1.shape = ", w1.shape
	w2 = ddn_helmholtz_fs(k, D_a, x1)
	print w1[0]
	print w2[0]
	assert w1[0] == w2[0]
	
	w2 = ddn_helmholtz_fs_j(k, D_a, x2, 0)
	print w2[0]
	print w2[-1]
	
def test2():
	a_radius = 0.1
	param1 = controlregion.make_annularsector_param_dict(0.11, 0.14, -np.pi/4, np.pi/4, 0, 16, 16, 4)
	param2 = controlregion.make_ball_param_dict(10, 64)
	D_c = controlregion.AnnularSector(param1)
	D_far = controlregion.Ball(param2)
	D_a = antenna.PolarArray(antenna.make_circle_s(a_radius), antenna.circle_ds, antenna.make_centers(np.cos, np.sin, 0, np.pi/2, 5), 32)
	k = 1
	
	u0 = boundaryfunction.PointSource(np.array([10000,0]))
	CP = ControlProblem(D_a, D_c, D_far, u0, k, 'delta')
	#print CP.A_near[:,0]
	#print CP.A_far[:,0]
	print CP.A_near.shape
	A = CP.A()
	print A[:,0]
	print A.shape
	B = np.dot(np.conjugate(A).T, A)
	print B.shape

	alpha = 1e-7
	Balpha = B + alpha*np.eye(B.shape[0])
	z = np.random.rand(B.shape[0])
	f = np.dot(B, z)
	x = np.linalg.solve(Balpha, f)
	print np.linalg.norm(np.dot(B, x) - f)
	print np.linalg.norm(x-z)
	print x.shape
	print x[0:20]
	print z[0:20]
	print np.dot(B,x)[0:20]
	print f[0:20]

def test3():
	a_radius = 0.1
	param1 = controlregion.make_annularsector_param_dict(0.11, 0.15, 3*np.pi/4, 5*np.pi/4, 0, 128, 128, 32)
	param2 = controlregion.make_ball_param_dict(10, 256)
	D_c = controlregion.AnnularSector(param1)
	#D_a = antenna.PolarArray(antenna.make_circle_s(a_radius), antenna.circle_ds, antenna.make_centers(np.cos, np.sin, 0, np.pi/2, 5), 32)
	D_a = antenna.Polar(antenna.make_circle_s(a_radius), antenna.circle_ds, np.array([0,0]), 256)
	k = 1
	D_far = controlregion.Ball(param2)
	u0 = boundaryfunction.PlaneWave(0)
	
	CP = ControlProblem(D_a, D_c, D_far, u0, k, 'fourier')
	A = CP.A()
	A_near = CP.A_near
	A_far = CP.A_far
	print(A_near[:,3])
	#print(D_c.get_boundary().T)
#	w = np.random.rand(D_a.n_points())
#	phi = CP.invert_coeff_vector(w)
#	print phi.shape
#	print len(phi)/D_a.n_points()
#	v = CP.control_trace(phi)
#	print v.shape
#	print v[0:D_a.n]
	
if __name__=="__main__":
	print "Running tests..."	
	#test1()
	#test2()
	test3()
