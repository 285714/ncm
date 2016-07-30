# stable periodic trajectories for ρ ∈ {350, 160, 100.5, 99.65, 28}
const σ=10, β=8/3
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


# plot the fixpoints and intersections
@noinline function cbPlotSolution(v)
	ρ = v[end]
	d = sqrt(β*(ρ-1))
	scatter3D([d;-d], [d;-d], ρ-1, edgecolors="none", c="r")

	ts,fs,ns = lorenzProjectionInternal(interps(unwrap(v)[1]),ρ,length(v)÷6*2)
	tmp = transpose(reduce(hcat, fs))
	scatter3D(tmp[:,1],tmp[:,2],tmp[:,3])
end


# takes 2π-periodic functions, n > Nyquist rate
# returns set of scalars to distinguish solutions
#TODO oversampling?
@noinline function lorenzProjectionInternal(ftmp,ρ,N)
	d = sqrt(β*(ρ-1))
	A,B = [d;d;ρ-1],[-d;-d;ρ-1]
	const p = 0.08872374069251765
	T = linspace(-pi+p,pi+p,2N)

	pA(t) = ftmp[1](t) + ftmp[2](t) + 2d #ftmp[1](t) - ftmp[2](t) #ftmp[3](t) - ρ + 1 ... interesting
	pB(t) = ftmp[1](t) + ftmp[2](t) - 2d
	t₀ = T[1]
	a₀ = pA(T[1])
	b₀ = pB(T[1])

	ts,fs,ns = [],[],[]
	for t in T[2:end]
		tmp = pA(t)
		if a₀ ≤ 0 < tmp
			push!(ts, bisection(pA, t₀, t))
			push!(fs, map(g->g(ts[end]), ftmp))
			push!(ns, norm(fs[end]-A))
		end
		a₀ = tmp

		tmp = pB(t)
		if b₀ ≤ 0 < tmp
			push!(ts, bisection(pB, t₀, t))
			push!(fs, map(g->g(ts[end]), ftmp))
			push!(ns, -norm(fs[end]-B))
		end
		b₀ = tmp

		t₀ = t
	end

	return ts, fs, sort(ns)
end

@noinline lorenzProjection(ftmp,ρ,N) = lorenzProjectionInternal(ftmp,ρ,N)[3]


#=
lorenzProjection(v...) = lorenzProjection2(v...) function lorenzProjection2(ftmp,ρ,N)
	T = linspace(-.1,2pi-.1,3N+2)
	plane(t) = ftmp[1](t) - ftmp[2](t)
	t₀,p₀ = T[1], plane(T[1])

	rtn = []
	for t in T[2:end-1]
		tmp = plane(t)
		p₀ ≤ 0 < tmp && push!(rtn, bisection(plane, t₀, t, ϵ=1e-4))
		t₀,p₀ = t,tmp
	end

	return map(x->norm(map(g->g(x), ftmp)), rtn)
end
=#

# N = 2m
# ρ = ℵ
# tmp = sqrt(β*(ρ-1))
# # fixpoints
# A,B = [tmp; tmp; ρ-1], [-tmp; -tmp; ρ-1]
#
# # vector from one fixpoint to the other, used as normal
# n = [ 2tmp; 2tmp; 0 ]
#
# T = linspace(-.1,2pi-.1,N+2)
# pA(t) = dot(n, f(T[1])-A)
# pB(t) = dot(n, f(T[1])-B)
# t₀ = T[1]
# a₀, b₀ = pA(t₀), pB(t₀)
#
# rtn = []
#
# I = T[2:end-1]
# state = start(I)
# t,state = next(I,state)
# tmp = pA(t)
# a₀ ≤ 0 < tmp && push!(rtn, bisection(pA, t₀, t, ϵ=1e-4))
# a₀ = tmp
#
# tmp = pB(t)
# b₀ ≤ 0 < tmp && push!(rtn, bisection(pB, t₀, t, ϵ=1e-4))
# b₀ = tmp
#
# t₀ = t
#
# rtn
#
# return map(f, rtn)

#=
function orthogonal_proj(zfront, zback)
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    return [
		1 0 0 0
		0 1 0 0
		0 0 a b
		0 0 0 zback]
end
proj3d.persp_transformation = orthogonal_proj
=#
