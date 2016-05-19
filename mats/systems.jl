σ,β,ρ = 10.0,8/3,100.0
lorenz(v::Vector) = [σ*(v[2]-v[1]); v[1]*(ρ-v[3])-v[2]; v[1]*v[2]-β*v[3]]
lorenz(t,v) = lorenz(v)
