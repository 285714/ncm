using PyPlot

include("lib/matsboUTIL.jl")
include("lib/matsboNWTN.jl")
include("lib/matsboRK.jl")
include("lib/matsboPRED.jl")

include("lib/ncmprojFINDINITIALDATA.jl")
include("lib/ncmprojPREPAREDATA.jl")
# include("lib/ncmprojGENERATEFUNCTIONS.jl")


include("systems.jl")
include("globals.jl")

include("callbacks.jl")


### Start ###
# TODO encapsulate, externalize, globalize...
#		t₀, y₀, (pred ~> PC...), H, J, rToC, cToR, ω, interpolate
# TODO ac border cases, transient plot
# TODO Samples prop to num periods
# NOTE Samples too small can lead to divergence

# System = lorenz
# TransientStepSize = .01
Samples = 256
c = 12.6

# H, J, rToC, cToR 					= generateFunctions(System, Dimensions, Epsilon)
Transient, SteadyState, P	= findCycle(System, .0, 10*rand(3), TransientIterations, TransientStepSize, SteadyStateIterations, SteadyStateStepSize)
y,ω												= prepareCycle(SteadyState, P, InterpPrec, Samples, fac=2)

# plots of approx steady-state trajectory, extracted cycle, resampled/shifted version
begin
	ioff(); figure(); ax = subplot(111,projection="3d"); hold(true)
	plot(Transient[:,1], Transient[:,2], Transient[:,3], color="y")
	plot(SteadyState[:,1], SteadyState[:,2], SteadyState[:,3], color="b")
	#plot(cyc[:,1], cyc[:,2], cyc[:,3], color="g")
	plot(y[:,1], y[:,2], y[:,3], color="r")
	ax["scatter"](y[1,1], y[1,2], y[1,3], color="r")
	show(); ion()
end

# TODO encapsulate...
Y = rfft(y, [1])
V₀ = [ real(Y[1,1]); real(Y[2:end,1]); imag(Y[2:end,1]); real(Y[1,2]); real(Y[2:end,2]); imag(Y[2:end,2]); real(Y[1,3]); real(Y[2:end,3]); imag(Y[2:end,3]); ω ]
V = matsboNWTN.newton(Hroessler, Jroessler, V₀, matsboPRED.predCount(16); init=initPlotResidual, callback=callbackPlotResidual)


Jtest = matsboNWTN.broyden(Hroessler, v -> matsboNWTN.forwardDifference(Hroessler, v))
Jtest = matsboNWTN.broyden(Hroessler, Jroessler)
V = matsboNWTN.newton(Hroessler, Jtest, V₀, matsboPRED.predCount(16); init=initPlotResidual, callback=callbackPlotResidual)
