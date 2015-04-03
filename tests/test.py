from pycloak import controlregion
import numpy as np

def main():
	param = {}
	param['r1'] = 0.11
	param['r2'] = 0.14
	param['theta1'] = 3*np.pi/4
	param['theta2'] = 5*np.pi/4
	param['doffset'] = 0.0
	param['narcin'] = 64
	param['narcout'] = 64
	param['ne'] = 32
	C1 = controlregion.AnnularSector(param)
	w = C1.get_weights()
	print(w)
	
if __name__ == "__main__":
	main()
