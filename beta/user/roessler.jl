# system files need to initialize global  Proj  variable

a,b = .1,.1
roessler(v) = [-v[2] - v[3]; v[1] + a*v[2]; b + v[3]*(v[1]-v[4])]
roessler(t,v) = roessler(v)

global Proj = Project()
