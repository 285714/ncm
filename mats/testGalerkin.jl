# lorenz,
# verzw
# interface
# defekt vs ortsraum
# startlösung homotopie

# TODO externalize predicates
# TODO externalize 'bifurcation indicator'
# TODO think of something for global ℵ
# TODO encapsulate, externalize, globalize...
# TODO Samples prop to num periods
# TODO step-length control feedback
# TODO PC criterion externalize
# TODO lib requirements
# TODO prevent autoscale
# NOTE Samples too small can lead to divergence

cd("/home/matsbo/Studies/ncm/mats")
push!(LOAD_PATH,"$(pwd())/lib")
using PyPlot

using matsboPRED
using matsboUTIL
using matsboNWTN
using matsboRK
using matsboINTERPOLATE
using matsboCONTINUATION

include("lib/ncmprojFINDINITIALDATA.jl")
include("lib/ncmprojPREPAREDATA.jl")
# include("lib/ncmprojGENERATEFUNCTIONS.jl")
include("lib/ncmprojCONTINUATION.jl")

include("globals.jl")
include("roessler.jl")
include("callbacks.jl")

### Start ###

ℵ = 4.0
InternSamples = 64

Transient, SteadyState, P	= findCycle(System, .0, 10*rand(3), TransientIterations, TransientStepSize, SteadyStateIterations, SteadyStateStepSize)
y,ω							= prepareCycle(SteadyState, P, fac=2)
V = [fromTraj(y); ω; ℵ]

H′(V) = H([V;ℵ]); J′(V) = J([V;ℵ])[:,1:end-1]
V = newton(H′, J′, V[1:end-1], predEps(1e-10)); V = [V;ℵ]

V = newton(H, J, V, predCount(8); init=initPlotResidual, callback=callbackPlotResidual)
V = newton(H, J, V, predCount(8); init=initPlot3D, callback=callbackPlot3D)
V = newton(H, J, V, predEps(1e-10))
norm(H(V))

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

# derivatives
# Jtest = broyden(H, v -> forwardDifference(H, v))
# Jtest = broyden(H, J)
# Jtest = v -> forwardDifference(H, v)
# V = newton(Hroessler, Jtest, V₀, predCount(16); init=initPlotResidual, callback=callbackPlotResidual)

global pl = [] #line array.. careful...
dir = true
W = PC(H, J, V, V -> 3.0 ≤ V[end] ≤ 20.0; init=initPC2, callback=callbackPC2, h₀=1.0, κ₀=0.5, δ₀=1.0, α₀=pi/180*2, dir=dir)

W = PC(V->H(V)+1e-2, J, W, V -> abs(V[end]-W[end]) ≤ 1; h₀=1.0, κ₀=0.5, δ₀=1.0, α₀=pi/180*2, dir=dir)
W = PC(H, J, W, V -> 3.0 ≤ V[end] ≤ 20.0; init=initPC2, callback=callbackPC2, h₀=1.0, κ₀=0.5, δ₀=1.0, α₀=pi/180*2, dir=dir)

W₁ = W
W₂ = W

figure(); subplot(111,projection="3d"); hold(true)
tmp = toTraj(W₁[1:end-2]); plot(tmp[:,1], tmp[:,2], tmp[:,3])
tmp = toTraj(W₂[1:end-2]); plot(tmp[:,1], tmp[:,2], tmp[:,3])

close()
