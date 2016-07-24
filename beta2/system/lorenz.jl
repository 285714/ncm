# stable periodic trajectories for ρ ∈ {350, 160, 100.5, 99.65, 28}
global σ=10, β=8/3
lorenz(t,v) = [
	σ*(v[2]-v[1])
	v[1]*(v[4]-v[3])-v[2]
	v[1]*v[2]-β*v[3]
	]

lorenz′(t,v) = [
	-σ σ 0 0
	(v[4]-v[3]) -1 -v[1] v[1]
	v[2] v[1] -β 0
	]
