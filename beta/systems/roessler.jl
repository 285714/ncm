# system files need to initialize global  Proj  variable

a,b,c = .1,.1,7
roessler(v) = [-v[2] - v[3]; v[1] + a*v[2]; b + v[3]*(v[1]-c)]
roessler(t,v) = roessler(v)

global Proj = Project(PC(), Galerkin(roessler, 3))
