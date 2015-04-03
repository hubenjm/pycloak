from pycloak import controlregion, antenna, dlpotential, regsolve
import pycloak.boundaryfunction as bdf	
import numpy as np

def main():

	k=10; alpha_i = 1.0; alpha_min = 1e-16; delta = 0.02; epsilon = 0.005; beta = 1.015; max_iter = 300;
	centers = np.array([0,0])
	x0 = np.array([10000,0])
	narcin = 256; narcout = 256; ne = 32; n = 256;
	r1 = 0.011; r2 = 0.014; theta1 = 3*np.pi/4; theta2 = 5*np.pi/4; doffset = 0;
	R = 10; nR = 256;
	
	D_a = antenna.PolarArray(antenna.make_circle_s(0.01), antenna.circle_ds, centers, n)
	
	param_a = controlregion.make_annularsector_param_dict(r1, r2, theta1, theta2, doffset, narcin, narcout, ne)
	param_R = controlregion.make_ball_param_dict(R, nR)
	D_c = controlregion.AnnularSector(param_a)
	B_R = controlregion.Ball(param_R)
	
	v1 = bdf.PointSource(x0)
	f1 = v1.trace(D_c, k)
	print "f1 has length", len(f1)
	
	f2 = np.zeros(nR, dtype=complex)
	print "f2 has length", len(f2)
	f = np.concatenate((f1,f2))
	
	CP = dlpotential.ControlProblem(D_a, D_c, B_R, v1, k, basis_type='delta')

	A = CP.A()
	print "A has shape", A.shape
	results = regsolve.tikhonov_solve(CP, alpha_i, alpha_min, delta, beta, max_iter, newton_tol = 1e-4, debug=False)
	
	print "near field relative error: {}".format(results['nfe'])
	print "far field L_inf error: {}".format(results['ffe'])
	print "norm of phi: {}".format(results['phinorm'])
	print "alpha: {}".format(results['alpha'])
	print "iterations: {}".format(results['iterations'])
	print "exit flag: {}".format(results['flag'])
		
if __name__ == "__main__":
	main()




