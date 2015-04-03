import numpy as np
from scipy.integrate import simps
from scipy.special import hankel1

class ControlRegion:
	
	def integrate(self, A):
		assert False
	def L2(self, f): #return square of L2 norm
		return self.integrate(f*np.conjugate(f))
	def L2_ave(self, f): #return square of L2 norm average
		return self.L2(f)/self.measure()	
	def get_boundary(self):
		assert False
	def n_points(self):
		assert False
	def supnorm(self, f):
		assert len(f) == self.n_points()
		return np.max(np.abs(f))
	def measure(self):
		assert False
		
	def planewavetrace(self, k, xi):
		xihat = np.array([np.cos(xi), np.sin(xi)])
		x = self.get_boundary()
		return np.exp(1j*k*(xihat[0]*x[0,:] + xihat[1]*x[1,:]))

	def pointsourcetrace(self, k, x0):
		assert x0.size == 2
		x = self.get_boundary()
		return (1j/4)*hankel1(0,k*np.linalg.norm(x - np.reshape(x0, (2,1)), axis=0))
		
		
class AnnularSector(ControlRegion):
	
	def __init__(self, param):
		assert param['theta2'] > param['theta1'] and param['r1'] < param['r2']
		self.param = {}
		self.param['r1'] = param['r1']
		self.param['r2'] = param['r2']
		self.param['theta1'] = param['theta1']
		self.param['theta2'] = param['theta2']
		self.param['doffset'] = param['doffset']
		self.param['narcin'] = param['narcin']
		self.param['narcout'] = param['narcout']
		self.param['ne'] = param['ne']
			
	def get_thetas(self, n, reverse=False):
		"""n is a positive integer representing the number of increments to use on an arc of the boundary. Returns a numpy array of shape (n,) whose values are the angles along the an arc of the control region. 
		"""
		dtheta = (self.param['theta2'] - self.param['theta1'])/n
		if reverse:
			return np.linspace(self.param['theta2'] - dtheta/2.0, dtheta/2.0 + self.param['theta1'], n)
		else:
			return np.linspace(dtheta/2.0 + self.param['theta1'], self.param['theta2'] - dtheta/2.0, n)
	
	def get_inner_points(self):
		thetas = self.get_thetas(self.param['narcin'])
		return np.array([self.param['doffset'] + self.param['r1']*np.cos(thetas), self.param['r1']*np.sin(thetas)])
	
	def get_outer_points(self):
		#should return points in clockwise order
		thetas = self.get_thetas(self.param['narcout'], reverse=True)
		return np.array([self.param['doffset'] + self.param['r2']*np.cos(thetas), self.param['r2']*np.sin(thetas)])
		
	def get_edgetop(self):
		ne = self.param['ne']
		ds_e = (self.param['r2'] - self.param['r1'])/ne
		r_vals = np.arange(self.param['r1'] + ds_e/2.0, self.param['r2'], ds_e)
		return np.reshape(np.array([self.param['doffset'] + np.cos(self.param['theta2'])*r_vals, np.sin(self.param['theta2'])*r_vals]), (2,len(r_vals)))
	
	def get_edgebottom(self, ne = None):
		#should return points in clockwise order
		if ne==None:
			ne = self.param['ne']
		ds_e = (self.param['r2'] - self.param['r1'])/ne
		r_vals = np.arange(self.param['r2'] - ds_e/2.0, self.param['r1'], -ds_e)
		return np.reshape(np.array([self.param['doffset'] + np.cos(self.param['theta1'])*r_vals, np.sin(self.param['theta1'])*r_vals]), (2, len(r_vals)))
		
	def get_boundary(self):
		return np.concatenate( (self.get_inner_points(), self.get_edgetop(), self.get_outer_points(), self.get_edgebottom()), axis=1)
		
	def get_weights(self):
		"""returns a numpy.array of shape (narcin + narcout + 2*ne, ) whose values are the quadrature weights to use when numerically integrating along the boundary of self
		"""
		narcin = self.param['narcin']
		narcout = self.param['narcout']
		ne = self.param['ne']
			
		theta_diff = self.param['theta2'] - self.param['theta1']
		dtheta_in = theta_diff/narcin
		dtheta_out = theta_diff/narcout
		ds_e = (self.param['r2'] - self.param['r1'])/ne
		return np.concatenate((np.repeat(self.param['r1']*dtheta_in, narcin), np.repeat(ds_e, ne), np.repeat(self.param['r2']*dtheta_out, narcout),  np.repeat(ds_e, ne)))
		
	def n_points(self):
		return self.param['narcin'] + self.param['narcout'] + 2*self.param['ne']
		
	def integrate(self, f):
		assert len(f) == self.n_points()
		return np.sum(f*self.get_weights())
	
	def supnorm(self, f):
		assert len(f) == self.n_points()
		return np.max(np.abs(f))
		
	def measure(self):
		dtheta = (self.param['theta2'] - self.param['theta1'])
		dr = self.param['r2'] - self.param['r1']
		return dtheta*(self.param['r1'] + self.param['r2']) + 2*dr
		
class Ball(ControlRegion):
	def __init__(self, param):
		self.param = {}
		self.param['R'] = param['R']
		self.param['n'] = param['n'] #number of sample points
		
	def integrate(self, f):
		assert len(f) == self.param['n']
		return sum(2*np.pi*self.param['R']/self.param['n']*f)
		
	def supnorm(self, f):
		assert len(f) == self.n_points()
		return np.max(np.abs(f))
	
	def get_thetas(self):
		return np.linspace(0, 2*np.pi - 2*np.pi/self.param['n'], self.param['n'])
		
	def get_boundary(self):
		thetas = self.get_thetas()
		return self.param['R']*np.array([np.cos(thetas), np.sin(thetas)])
	
	def n_points(self):
		return self.param['n']
		
	def measure(self):
		return 2*np.pi*self.param['R']
		
class Rectangle(ControlRegion):
	def __init__(self, param):
		assert param['xll'].shape == (2,)
		self.param = {}
		self.param['xll'] = param['xll'] #coordinates of lower left corner of rectangle
		self.param['height'] = param['height']
		self.param['width'] = param['width']
		self.param['n'] = param['n']
		
	def n_points(self):
		return self.param['n']
		
	def integrate(self, f):
		assert len(f) == self.n_points()
		return np.sum(self.ds()*f)
		
	def measure(self):
		return 2*self.param['height'] + 2*self.param['width']
		
	def supnorm(self, f):
		assert len(f) == self.n_points()
		return np.max(np.abs(f))
	
	def ds(self):
		return self.measure()/self.n_points()
		
	def get_boundary(self):
		"""Returns a 2 by n matrix whose columns are evenly spaced points along the antenna with respect to arclength.
		"""
		xll = self.param['xll']
		L = np.linspace(0, self.measure() - self.ds(), self.n_points())
		j1 = len(L[L<self.param['height']])
		j2 = len(L[L<self.param['height'] + self.param['width']])
		j3 = len(L[L<2*self.param['height'] + self.param['width']])
		j4 = self.n_points()

		a1 = xll[0]*np.ones(j1)
		a2 = (xll[0] + self.param['width'])*np.ones(j3-j2)
		b1 = (xll[1] + self.param['height'])*np.ones(j2-j1)
		b2 = (xll[1] - self.param['height'])*np.ones(j4-j3)

		X = np.hstack([a1, L[j1:j2]-self.param['height'] - self.param['width']/2.0, a2, -L[j3:]+2*self.param['height']+1.5*self.param['width']])
		Y = np.hstack([L[:j1]-self.param['height']/2.0, b1, -L[j2:j3]+1.5*self.param['height'] + self.param['width'], b2])
		return np.vstack([X,Y])

def make_annularsector_param_dict(r1, r2, theta1, theta2, doffset, narcin, narcout=None, ne=None):
	param = {}
	if not narcout:
		dtheta = theta2-theta1
		narcout = int(np.ceil(r2*narcin/r1))
	if not ne:
		dtheta = theta2-theta1
		dr = r2-r1
		ne = int(np.ceil(dr*narcin/(r1*dtheta)))
	param['r1'] = r1; param['r2'] = r2; param['theta1'] = theta1; param['theta2'] = theta2
	param['doffset'] = doffset; param['narcin'] = narcin; param['narcout'] = narcout; param['ne'] = ne
	return param

def make_rectangle_param_dict(height, width, xll, n):
	param = {}
	param['height'] = height
	param['width'] = width
	param['xll'] = xll
	param['n'] = n
	return param
	
def make_ball_param_dict(R, n):
	param = {}
	param['R'] = R
	param['n'] = n
	return param
	
def test1():
	param = {}
	param['r1'] = 0.11; param['r2'] = 0.14; param['theta1'] = -np.pi/4; param['theta2'] = np.pi/4
	param['doffset'] = 0.0; param['narcin'] = 64; param['narcout'] = 64; param['ne'] = 32
	C1 = AnnularSector(param)
	w = C1.get_weights()
	print w, w.shape
	print C1.param['theta1']
	print C1.param['narcin'], C1.param['narcout'], C1.param['ne']
	
	
	paramrect = {}
	paramrect['height'] = 0.5; paramrect['width'] = 0.1; paramrect['xll'] = np.array([0.11, -0.25])
	paramrect['n'] = 256
	C2 = Rectangle(paramrect)
	
	paramball = {}
	paramball['R'] = 100; paramball['n'] = 512
	C3 = Ball(paramball)
	
	N = C2.n_points() + C3.n_points()
	f = np.random.rand(N)
	g = np.random.rand(N)
	print C1.L2(f[:C1.n_points()])
	
if __name__=="__main__":
	print "Running tests..."	
	test1()
	
