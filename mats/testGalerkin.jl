using PyPlot

include("lib/matsboUTIL.jl")
include("lib/matsboNWTN.jl")
include("lib/matsboRK.jl")
include("lib/matsboPRED.jl")

include("lib/ncmprojFINDINITIALDATA.jl")
include("lib/ncmprojPREPAREDATA.jl")
include("lib/ncmprojGENERATEFUNCTIONS.jl")

include("systems.jl")
include("globals.jl")


# callback for interactive plotting
global ax1, ax2, distance
function initPlot()
		global ax1, ax2, distance
		distance = Float64[]
    ion(); figure()
		global ax1 = subplot2grid((4,1), (0, 0), rowspan=3, projection="3d"); plot([.0],[.0],[.0])
		global ax2 = subplot2grid((4,1), (3, 0))
end

function callbackPlot(V, H, J)
	global ax1, ax2, distance

  local C,ω,N
  C,ω = rToC(V)
	N = 2*(size(C,1)-1)
  sca(ax1); ax1["lines"][end]["set_color"]("gray"); hold(true); plot(irfft(C[:,1],N), irfft(C[:,2],N), irfft(C[:,3],N), color="r")

	push!(distance, sumabs(H))
	sca(ax2); hold(false); plot(distance)

	draw(); sleep(.05)
end



### Start ###
# TODO encapsulate, externalize, globalize...
#		t₀, y₀, (pred ~> PC...), H, J, rToC, cToR, ω, interpolate
# TODO ac border cases, transient plot

H, J, rToC, cToR 		= generateFunctions(System, Dimensions, Epsilon)
Trace, P	 					= findCycle(System, .0, 10*rand(3), TransientIterations, TransientStepSize, SteadyStateIterations, SteadyStateStepSize)
y,ω									= prepareCycle(Trace, P, Samples, fac=2)

# plots of approx steady-state trajectory, extracted cycle, resampled/shifted version
ioff(); figure(); ax = subplot(111,projection="3d"); hold(true)
plot(Trace[:,1], Trace[:,2], Trace[:,3], color="b")
#plot(cyc[:,1], cyc[:,2], cyc[:,3], color="g")
plot(y[:,1], y[:,2], y[:,3], color="r")
ax["scatter"](y[1,1], y[1,2], y[1,3], color="r")
show()

# TODO encapsulate...
Y = rfft(y, [1])
V₀ = cToR(Any[matsboUTIL.columns(Y); ω])
V = matsboNWTN.newton(H, matsboNWTN.broyden(H, J), V₀, matsboPRED.predCount(16); init=initPlot, callback=callbackPlot)
# V = matsboNWTN.newton(H, J, V₀, matsboPRED.predCount(8); init=initPlot, callback=callbackPlot)

# optimize
# @time V₁ = matsboNWTN.newton(H, J, V₀, matsboPRED.predCount(16); init=initPlot, callback=callbackPlot)
# @time V₂ = matsboNWTN.newton(H, matsboNWTN.broyden(H, J), V₀, matsboPRED.predCount(32); init=initPlot, callback=callbackPlot)

norm(H(V))
norm(H(V₀))
# norm(H(V₁))
# norm(H(V₂))

# sum(real(Y[:,1]))
