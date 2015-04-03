import numpy as np
from pycloak import controlregion

n = 200 #number of antenna sample pts
m = 300 #number of control region sample pts

A = np.random.rand(m,n) + 1j*np.random.rand(m,n)
xexact = np.random.rand(n) + 1j*np.random.rand(n)
yexact = np.dot(A,xexact)

w = np.dot(np.conjugate(A.T), yexact)
B = np.dot(np.conjugate(A.T), A) + 0.01*np.eye(n)
x1 = np.linalg.solve(B,w)
print np.linalg.norm(x1-xexact)
print np.linalg.norm(np.dot(A,x1) - yexact)
print xexact[0:10]
print x1[0:10]

A = np.random.rand(m,n)
xexact = np.random.rand(n)
yexact = np.dot(A,xexact)

w = np.dot(np.conjugate(A.T), yexact)
B = np.dot(np.conjugate(A.T), A) + 0.01*np.eye(n)
x1 = np.linalg.solve(B,w)
print np.linalg.norm(x1-xexact)
print np.linalg.norm(np.dot(A,x1) - yexact)
print xexact[0:10]
print x1[0:10]

#the number of sample points on ControlRegion should be
#AT LEAST as big as the number of sample points on AntennaArray



