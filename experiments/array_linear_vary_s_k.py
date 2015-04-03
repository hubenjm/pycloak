from pycloak import controlregion, antenna, boundaryfunction, dlpotential, regsolve
import numpy as np
from plot import plot_s_k

#varies the spacing of the antennas in a linear array, oriented vertically

def vary_s_k(outfile, CP, s_min, s_max, N_s, k_min, k_max, N_k, alpha_0, alpha_min, delta, beta, max_iter, newton_tol = 1e-4, debug=False):
	"""
	"""
	#this is only currently intended for use with rectangular control regions
	assert s_min <= s_max
	assert k_min <= k_max and 0 < k_min
	
	s = np.linspace(s_min, s_max, N_s)
	k = np.linspace(k_min, k_max, N_k)
	
	#compute radius of antennae (constant in this experiment)
	a = np.max(CP.D_a.eval_s())
	
	#get number of sample points on each antennae in array
	n = CP.D_a.n
	
	#initialize data arrays
	NFE = np.zeros((N_s, N_k))
	FFE = np.zeros((N_s, N_k))
	F = np.zeros((N_s, N_k))
	PHINORM = np.zeros((N_s, N_k))
	FLAG = np.zeros((N_s, N_k))
	ITER = np.zeros((N_s, N_k))
	ALPHA = np.zeros((N_s, N_k))
	
	for i in range(N_s):
	#adjust number of antennae so that width of array equals width of rectangular control region D_c
		CP.set_D_a(compute_D_a(CP.D_c, s[i], a, n))
		for j in range(N_k):
			CP.set_k(k[j])
			results = regsolve.tikhonov_solve(CP, alpha_0, alpha_min, delta, beta, max_iter, newton_tol, debug)
			NFE[i,j] = results['nfe']
			FFE[i,j] = results['ffe']
			F[i,j] = results['F']
			PHINORM[i,j] = results['phinorm']
			FLAG[i,j] = results['flag'][0]
			ITER[i,j] = results['iterations']
			ALPHA[i,j] = results['alpha']
			print("s = {}, n_antenna = {}, k = {}, nfe = {}, ffe = {}, phinorm = {}, alpha = {}, flag = {}").format(s[i], CP.D_a.n_antennas, k[j], NFE[i,j], FFE[i,j], PHINORM[i,j], ALPHA[i,j], FLAG[i,j])
	
	#build 2D arrays for k and s
	K = np.tile(k, (N_s, 1))
	S = np.tile(s, (N_k, 1)).T

	np.savez(outfile, DELTA = delta, S=S, K=K, NFE=NFE, FFE=FFE, F=F, PHINORM=PHINORM, FLAG=FLAG, ITER=ITER, ALPHA=ALPHA)
	

def compute_D_a(D_c, s, a, n):
	"""
	D_c = control region object which we want to align the antenna array with
	s = spacing between antennae
	a = radius of each antenna
	n = number of sample points on each antenna
	"""
	assert isinstance(D_c, controlregion.Rectangle)
	mid = D_c.param['height']/2
	N_a = int(2*np.max([0, np.floor((mid-a-(s/2))/(a+s))])) + 1
	
	centers = np.array([np.zeros(N_a), np.linspace(-mid, mid, N_a)])
	D_a = antenna.PolarArray( antenna.make_circle_s(a), antenna.circle_ds, centers, n)
	return D_a

def main():
	
	outfile = "array_linear_vary_s_k_data"
	a = 0.01
	s_min = 0.015
	s_max = 0.02
	N_s = 20
	k_min = 1
	k_max = 20
	N_k = 20
	alpha_0 = 1.0
	alpha_min = 1e-16
	delta = 0.02
	beta = 1.015
	max_iter = 300

	height = 0.2
	width = 0.01
	xll = np.array([-0.021, -0.1])
	n = 32
	
	n_c = 128
	n_R = 128
	R = 10.0
	x_0 = np.array([1000,0])
	
	param = controlregion.make_rectangle_param_dict(height, width, xll, n)
	param_R = controlregion.make_ball_param_dict(R, n_R)
	D_c = controlregion.Rectangle(param)
	D_a = antenna.Polar(antenna.make_circle_s(a), antenna.circle_ds, np.array([0,0]), n)
	D_far = controlregion.Ball(param_R)
	u0 = boundaryfunction.PointSource(x_0)
	k = k_min
	
	CP = dlpotential.ControlProblem(D_a, D_c, D_far, u0, k, 'delta')
	vary_s_k(outfile, CP, s_min, s_max, N_s, k_min, k_max, N_k, alpha_0, alpha_min, delta, beta, max_iter)
	

if __name__ == "__main__":
	main()

