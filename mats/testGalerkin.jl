# TODO externalize predicates
# TODO externalize 'bifurcation indicator'
# TODO think of something for global ℵ
# TODO encapsulate, externalize, globalize...
# TODO Samples prop to num periods
# NOTE Samples too small can lead to divergence

cd("/home/matsbo/Studies/ncm/mats")
push!(LOAD_PATH,"$(pwd())/lib")

using PyPlot
using matsboPRED

include("lib/matsboUTIL.jl")
include("lib/matsboNWTN.jl")
include("lib/matsboRK.jl")
#include("lib/matsboPRED.jl")

include("lib/ncmprojFINDINITIALDATA.jl")
include("lib/ncmprojPREPAREDATA.jl")
# include("lib/ncmprojGENERATEFUNCTIONS.jl")
include("lib/ncmprojCONTINUATION.jl")

include("globals.jl")
include("roessler.jl")

include("callbacks.jl")

### Start ###

ℵ = 12.0
InternSamples = 128

Transient, SteadyState, P	= findCycle(System, .0, 10*rand(3), TransientIterations, TransientStepSize, SteadyStateIterations, SteadyStateStepSize)
y,ω							= prepareCycle(SteadyState, P, fac=1)
V = [fromTraj(y); ω; ℵ]


# H′(V) = H([V,ℵ]); J′(V) = J([V,ℵ])[:,1:end-1]
# V = matsboNWTN.newton(H′, J′, V[1:end-1], matsboPRED.predCount(32); init=initPlotResidual, callback=callbackPlotResidual); V = [V;ℵ]

V = matsboNWTN.newton(H, J, V, predCount(48)∧predEps(1e-6))
V = matsboNWTN.newton(H, J, V, matsboPRED.predCount(16); init=initPlotResidual, callback=callbackPlotResidual)
V = matsboNWTN.newton(H, J, V, matsboPRED.predCount(16); init=initPlot3D, callback=callbackPlot3D)
V = matsboNWTN.newton(H, J, V, predEps(1e-2))
norm(H(V))

tmp = V
function cb(v,h,j)
	global tmp
	println(norm(v-tmp))
	tmp = v
end
V = matsboNWTN.newton(H, J, V, matsboPRED.predCount(32); callback=cb)

B,W,ℵ′ = incrementalLoading(V, .1, 3.0)
B,W,ℵ′ = incrementalLoading(V, .1, 18.0)
# B,W,ℵ′ = incrementalLoading(V, .05, ℵ-1; init=initIncLoad, callback=callbackIncLoad)

# bifurcation diagrams
tmpc = map(first, B); c="b"
for i in 2:maximum(map(length, B))
	tmp = map(x->i≤length(x)?x[i]:.0, B)
	hold(true)
	plot(tmpc,tmp,color=c)
end

gca()["lines"][end]["remove"]()

# T, SS, ...
begin
	local y = toTraj(V[1:end-2])
	figure(); ax = subplot(111,projection="3d"); hold(true)
	plot(Transient[:,1], Transient[:,2], Transient[:,3], color="y")
	plot(SteadyState[:,1], SteadyState[:,2], SteadyState[:,3], color="b")
	#plot(cyc[:,1], cyc[:,2], cyc[:,3], color="g")
	plot(y[:,1], y[:,2], y[:,3], color="r")
	ax["scatter"](y[1,1], y[1,2], y[1,3], color="r")
end

# Traj before/after cont
figure(); subplot(111,projection="3d"); hold(true)
tmp = toTraj(V[1:end-2]); plot(tmp[:,1], tmp[:,2], tmp[:,3])
tmp = toTraj(W[1:end-2]); plot(tmp[:,1], tmp[:,2], tmp[:,3])

# derivatives
# Jtest = matsboNWTN.broyden(H, v -> matsboNWTN.forwardDifference(H, v))
# Jtest = matsboNWTN.broyden(H, J)
# Jtest = v -> matsboNWTN.forwardDifference(H, v)
# V = matsboNWTN.newton(Hroessler, Jtest, V₀, matsboPRED.predCount(16); init=initPlotResidual, callback=callbackPlotResidual)

# System = (t,v) -> roessler(t,[v; ℵ])
# H = v -> Hroessler([v; ℵ])
# J = v -> Jroessler([v; ℵ])[:,1:end-1]

function initPC()
	global P
	hold(true)
	ion()
	P==Void && (P = scatter([], [], s=1, color="k"))
	# show()
	# sleep(.01)
end

function callbackPC(V)
	global P
	local v = toTraj(V[1:end-2])
	local ofs = P["get_offsets"]()

	for i in 1:size(v,1)-1
		v[i,1] ≤ 0 && v[i+1,1] > 0 && (ofs = vcat(ofs, [V[end] norm(v[i,:])]))
	end
	v[end,1] ≤ 0 && v[1,1] > 0 && (ofs = vcat(ofs, [V[end] norm(v[end,:])]))

	P["set_offsets"](ofs)
	gca()["set_xlim"](minimum(ofs[:,1])-.5, maximum(ofs[:,1])+.5)
	gca()["set_ylim"](minimum(ofs[:,2])-.5, maximum(ofs[:,2])+.5)
	# draw()
	# sleep(.01)

	return Void
end


P = Void
ℵ = V[end]
W = PC(H, J, V, V -> 3.0 ≤ V[end] ≤ 14.1; init=initPC, callback=callbackPC, h₀=1.0, κ₀=.5, δ₀=.1, dir=true)
W[end]

V=W
close()
