# stable periodic trajectories for ρ ∈ {350, 160, 100.5, 9.65}
include("../lib/ncmprojGALERKINHELPER.jl")

global σ=10, β=8/3
lorenz(t,v) = [
	σ*(v[2]-v[1])
	v[1]*(v[4]-v[3])-v[2]
	v[1]*v[2]-β*v[3]
	]
lorenz(v) = lorenz(0, v)

lorenz′(t,v) = [
	-σ σ 0 0
	(v[4]-v[3]) -1 -v[1] v[1]
	v[2] v[1] -β 0
	]
lorenz′(v) = lorenz′(0, v)

f, f′, H, J = lorenz, lorenz′, fToH(lorenz), f′ToJ(lorenz′)
