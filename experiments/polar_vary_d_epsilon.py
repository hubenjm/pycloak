from pycloak import controlregion, antenna, boundaryfunction, dlpotential, regsolve
import numpy as np

#varies the spacing of the antennas in a linear array, oriented vertically

def vary_d_epsilon(outfile, a, doffset, dr, theta1, theta2, R, k, n_a, n_c, n_R, x_0, d_min, d_max, N_d, epsilon_max, N_epsilon, alpha_0, alpha_min, delta, beta, max_iter, newton_tol = 1e-4, debug=False):
	"""
	"""
	#string format for output
	HEADERFORMAT = "   d \t \t epsilon \t    nfe \t    ffe \t   phinorm \t  stability \t   alpha \t flag"
	FORMAT = "{:.3e} \t {:.4f} \t {:.4e} \t {:.4e} \t {:.4e} \t {:.4e} \t {:.4e} \t  {:.0f}"
	
	#this is only currently intended for use with rectangular control regions
	assert 0 < d_min and d_min <= d_max
	assert 0 < epsilon_max
	
	d = np.linspace(d_min, d_max, N_d)
	epsilon = np.linspace(0, epsilon_max, N_epsilon)
	r1 = a + d - doffset
	r2 = r1 + dr
	param_R = controlregion.make_ball_param_dict(R, n_R) #doesn't change
	D_a = antenna.Polar(antenna.make_circle_s(a), antenna.circle_ds, np.array([0,0]), n_a) #doesn't change
	D_far = controlregion.Ball(param_R)
	u0 = boundaryfunction.PointSource(x_0)
	
	#initialize data arrays
	NFE = np.zeros((N_d, N_epsilon))
	FFE = np.zeros((N_d, N_epsilon))
	F = np.zeros((N_d, N_epsilon))
	PHINORM = np.zeros((N_d, N_epsilon))
	STABILITY = np.zeros((N_d, N_epsilon)) # ||phi_epsilon - phi_0||/||phi_0||
	FLAG = np.zeros((N_d, N_epsilon))
	ITER = np.zeros((N_d, N_epsilon))
	ALPHA = np.zeros((N_d, N_epsilon))
	
	#output header
	print(HEADERFORMAT)
	for i in range(N_d):
		param = controlregion.make_annularsector_param_dict(r1[i], r2[i], theta1, theta2, doffset, n_c)
		D_c = controlregion.AnnularSector(param)
		CP = dlpotential.ControlProblem(D_a, D_c, D_far, u0, k, 'fourier')
		
		#compute noiseless solution
		results = regsolve.tikhonov_solve(CP, alpha_0, alpha_min, delta, beta, max_iter, None, newton_tol, debug)
		NFE[i,0] = results['nfe']
		FFE[i,0] = results['ffe']
		F[i,0] = results['F']
		PHINORM[i,0] = results['phinorm']
		FLAG[i,0] = results['flag'][0]
		ITER[i,0] = results['iterations']
		ALPHA[i,0] = results['alpha']
		
		phi0 = CP.invert_coeff_vector(results['phi']) #solution for unperturbed data
		
		print('\n'+FORMAT).format(d[i], epsilon[0], NFE[i,0], FFE[i,0], PHINORM[i,0], 0, ALPHA[i,0], FLAG[i,0])
		
		#compute f1norm
		f1 = u0.trace(D_c, k)
		f1norm = np.sqrt(D_c.L2(f1))
		
		#compute noise vector
		nu = np.random.rand(D_c.n_points()) + 1j*np.random.rand(D_c.n_points())
		nu = nu/np.sqrt(D_c.L2(nu)) #normalize to unit norm
		noise = np.concatenate((nu, np.zeros(n_R)))
		
		for j in range(1, N_epsilon):
			#add noise epsilon[j]*noise
			results = regsolve.tikhonov_solve(CP, alpha_0, alpha_min, delta, beta, max_iter, epsilon[j]*f1norm*noise, newton_tol, debug)
			NFE[i,j] = results['nfe']
			FFE[i,j] = results['ffe']
			F[i,j] = results['F']
			PHINORM[i,j] = results['phinorm']
			FLAG[i,j] = results['flag'][0]
			ITER[i,j] = results['iterations']
			ALPHA[i,j] = results['alpha']
			
			phi1 = CP.invert_coeff_vector(results['phi'])
			STABILITY[i,j] = np.real(D_a.norm(phi1-phi0)/D_a.norm(phi0))
			
			print(FORMAT).format(d[i], epsilon[j], NFE[i,j], FFE[i,j], PHINORM[i,j], STABILITY[i,j], ALPHA[i,j], FLAG[i,j])
	
	#build 2D arrays for d and epsilon
	EPSILON = np.tile(epsilon, (N_d, 1))
	D = np.tile(d, (N_epsilon, 1)).T

	np.savez(outfile, DELTA = delta, D=D, EPSILON=EPSILON, NFE=NFE, FFE=FFE, F=F, PHINORM=PHINORM, STABILITY=STABILITY, FLAG=FLAG, ITER=ITER, ALPHA=ALPHA)


def main():
	
	outfile = "polar_vary_d_epsilon_data"
	a = 0.01
	dr = 0.004
	theta1 = 3*np.pi/4
	theta2 = 5*np.pi/4
	doffset = 0.0
	R = 10.0
	k = 1
	d_min = 0.001
	d_max = 0.02
	N_d = 30
	alpha_0 = 1.0
	alpha_min = 1e-16
	epsilon_max = 0.01
	N_epsilon = 20
	delta = 0.02
	beta = 1.015
	max_iter = 100

	n_a = 64	
	n_c = 64
	n_R = 64
	x_0 = np.array([10000,0])
	
	vary_d_epsilon(outfile, a, doffset, dr, theta1, theta2, R, k, n_a, n_c, n_R, x_0, d_min, d_max, N_d, epsilon_max, N_epsilon, alpha_0, alpha_min, delta, beta, max_iter)
	

if __name__ == "__main__":
	main()

