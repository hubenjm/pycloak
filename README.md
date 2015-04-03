PyCloak
====
Install
-------
Check out the master branch from the [pycloak][pycloakrepo].  Then,
    
    python setup.py install

[pycloakrepo]: http://www.github.com/hubenjm/pycloak

Description
-----------
This package is designed to solve the a particular active control problem in two dimensions for the homogeneous Helmholtz equation. Mathematical details are presented in the included pdf file, *./mathdoc/readme-overview.pdf*. This software was written initially in order to yield numerical results for a recent research paper. As such, it includes a variety of modules and routines designed to set up the general framework of the problem, allowing for manipulation of many different domain parameters and coefficients. Once the problem is set up properly using the ```ControlProblem``` class, one can use tikhonov_solve from regsolve.py in order to find a given solution.

API
-------

pycloak
1. antenna.py
	- Antenna - Generic class for an antenna object
	- Rectangle - Rectangular shaped antenna object
	- PolarArray - Multiple antenna array object, where each antenna has the same parametrization in polar coordinates with a different center.
	- Polar - A single antenna defined in polar coordinates and with center at the origin
	
2. controlregion.py
	- 
3. dlpotential.py
4. boundaryfunction.py
5. regsolve.py

### Deprecated / Incomplete ###

1. guicloaker.py
2. visuals.py


