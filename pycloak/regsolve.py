import numpy as np
import dlpotential
import antenna
import controlregion
import boundaryfunction
	
def tikhonov_solve(CP, alpha_0, alpha_min, delta, beta, max_iter, noise = None, newton_tol = 1e-4, debug=False):
	"""
	CP is a ControlProblem
`	alpha_0 is the initial regularization parameter
	alpha_min is the smallest allowed value of alpha (e.g. 1e-20)
	delta is the discrepancy to use when trying to solve ||Aphi - u0||^2 - delta^2 = 0
	max_iter is a positive integer bound on the number of iterations for the Newton Method
	beta is the multiplicative stepsize for the initial line search
	noise is a vector of noise (default is NoneType)
	"""
	
	alpha = float(alpha_0)
	A = CP.A()
	B_near = CP.B_near()
	B_far = CP.B_far()
	B = CP.B()
	A_near = CP.compute_A_near()
	A_far = CP.compute_A_far()
	f_near = CP.u0.trace(CP.D_c, CP.k)
	f_near_normsq = CP.D_c.L2(f_near)
	f_far = np.zeros(CP.D_far.n_points())
	f = np.concatenate((f_near, f_far), axis=0)
	
	#add noise if supplied
	if noise is not None:
		assert len(noise) == len(f)
		f = f + noise
		
	N = B.shape[0]
	g = np.dot(np.conjugate(A_near).T, f_near)
	h = np.linalg.solve(B + alpha*np.eye(N), g)
	F = CP.F(h, f, delta)
	if debug: print("initial error = {0}").format(np.real(F))
	while F > 0 and alpha > alpha_min:
		alpha /= beta
		h = np.linalg.solve(B + alpha*np.eye(N), g)
		F = CP.F(h, f, delta)
	
	errorsq = F + delta**2
	if debug:
		print("finished line search.")
		print("(nfe, ffe) = {0}, {1}").format(CP.nfe(h,f_near), CP.ffe(h,f_far))
		print("starting Newton's Method...")
	
	j = 0
	if F <= 0:
		while np.sqrt(np.abs(F)) > delta*newton_tol and j < max_iter:
			j += 1
			dh = -np.linalg.solve(B + alpha*np.eye(N), h)
			dF = -2.0/f_near_normsq*alpha*np.real(CP.D_a.integrate(np.conjugate(dh)*h))\
			+ 2*(1.0/CP.D_far.measure() - 1.0/f_near_normsq)\
			*np.real(CP.D_a.integrate(np.conjugate(dh)*(np.dot(B_far, h) - np.dot(np.conjugate(A_far).T, f_far))))
		
			temp = alpha - (F/dF)
			if temp > 0:
				exit_flag = (0, "Tikhonov regularization via Newton's method completed with no issues.")
				h_new = np.linalg.solve(B + temp*np.eye(N),g)
				F_new = CP.F(h_new, f, delta)
				errornewsq = F_new + delta**2
				if np.abs(F_new) <= np.abs(F):
					h = h_new
					errorsq = errornewsq
					alpha = temp;
				else:
					break
					exit_flag = (1, "Tikhonov regularization terminated during Newton's method due to non-decreasing F.")

			else:
				exit_flag = (2, "Tikhonov regularization terminated during Newton's method due to negative alpha after descent step.")
				break
				
			F = errorsq - delta**2
	else:
		exit_flag = (3, "Tikhonov regularization terminated before Newton's method due to no feasible value of alpha existing which is larger than minimum allowable alpha specified.")
	
	phi = CP.invert_coeff_vector(h)	
	results = {}
	results['phi'] = phi
	results['nfe'] = np.real(CP.nfe(h,f_near))
	results['ffe'] = np.real(CP.ffe(h,f_far))
	results['phinorm'] = np.real(CP.D_a.norm(phi))
	results['alpha'] = np.real(alpha)
	results['iterations'] = j
	results['flag'] = exit_flag
	results['F'] = F
	return results
	
def test1():
	alpha_0 = 1; alpha_min = 1e-20; k=1; a_radius = 0.01
	delta = 0.005; beta = 1.015; max_iter = 300;
	debug = True; newton_tol=1e-4;
	
	centers = np.array([[0,0,0,0,0],[-0.6, -0.3, 0.0, 0.3, 0.6]])
	
	D_a = antenna.PolarArray(antenna.make_circle_s(a_radius), antenna.circle_ds, centers, 64)
	D_a2 = antenna.Polar(antenna.make_circle_s(a_radius), antenna.circle_ds, np.array([0,0]), 256)
	param1 = controlregion.make_annularsector_param_dict(0.011, 0.016, 3*np.pi/4, 5*np.pi/4, 0, 128, 128, 32)
	param2 = controlregion.make_ball_param_dict(100, 256)
	D_c = controlregion.AnnularSector(param1)
	D_far = controlregion.Ball(param2)
	
	u0 = boundaryfunction.PlaneWave(np.pi/2.0)
	u1 = boundaryfunction.PointSource(np.array([10000,0]))
	
	CP1 = dlpotential.ControlProblem(D_a, D_c, D_far, u1, k, 'delta')
	results = tikhonov_solve(CP1, alpha_0, alpha_min, delta, beta, max_iter, None, newton_tol, debug)
	for j in results:
		if j != 'phi':
			print (j + ": {}").format(results[j])
			
if __name__=="__main__":
	print "Running tests..."	
	test1()	
	
	
	
	
	
	
