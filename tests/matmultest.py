import numpy as np

n = 4097
A = np.random.rand(n,n)
D = np.random.rand(n)
B = np.zeros((n,n))

def test1():
	for i in range(n):
		for j in range(n):
			B[i,j] = D[i]*A[i,j]
	
	print B[0,0]
	
if __name__=="__main__":
	print "Running tests..."	
	test1()
