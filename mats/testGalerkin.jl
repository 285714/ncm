# TODO externalize predicates
# TODO externalize 'bifurcation indicator'
# TODO think of something for global ℵ

using PyPlot

include("lib/matsboUTIL.jl")
include("lib/matsboNWTN.jl")
include("lib/matsboRK.jl")
include("lib/matsboPRED.jl")

include("lib/ncmprojFINDINITIALDATA.jl")
include("lib/ncmprojPREPAREDATA.jl")
# include("lib/ncmprojGENERATEFUNCTIONS.jl")

include("globals.jl")
include("roessler.jl")

include("callbacks.jl")

ℵ = 12.0
InternSamples = 128

### Start ###
# TODO encapsulate, externalize, globalize...
#		t₀, y₀, (pred ~> PC...), H, J, rToC, cToR, ω, interpolate
# TODO Samples prop to num periods
# NOTE Samples too small can lead to divergence


function toTraj(V)
	global ExternSamples
	local m = (length(V)-Dimensions)÷(2*Dimensions)
	local tmp = reshape(V, 2m+1, Dimensions)
	local P = map(1:Dimensions) do i
		x -> matsboUTIL.interpolateTrigonometric(tmp[1,i], 2tmp[2:2+m-1,i], -2tmp[2+m:end,i])(x) / (2m+1)
	end
	local x = linspace(.0, 2pi, ExternSamples)
	return reduce(hcat, map(f->f(x), P))
end

function fromTraj(v)
	global InterpPrec, InternSamples
	local y = mapslices(v->matsboUTIL.interpolate(v, InterpPrec, linspace(1.0, size(v,1), InternSamples)), v, [1])
	local Y = rfft(y, [1])
	return vec(mapslices(y -> [real(y); imag(y[2:end])], Y, [1]))
end

function incrementalLoading(V, Δℵ, ℶ; init=Void, callback=Void)
	global ℵ
	local B = Array{Float64}[]
	local W = deepcopy(V)
	local ℵ′ = ℵ

	init≠Void && init()
	callback≠Void && callback(ℵ, W, Void)

	for ℵ′′ in ℵ:(ℵ≤ℶ?Δℵ:-Δℵ):ℶ 	# =)
		ℵ = ℵ′′
		W = matsboNWTN.newton(H, J, W, matsboPRED.predCount(8))
		local w = toTraj(W[1:end-1])
		callback≠Void&& callback(ℵ, W, w)
		local maxs = Float64[]
		for i in 1:size(w,1)-1
			w[i,1] ≤ 0 && w[i+1,1] > 0 && push!(maxs, norm(w[i,:]))
		end
		w[end,1] ≤ 0 && w[1,1] > 0 && push!(maxs, norm(w[end,:]))
		(ℵ′<ℶ ? push! : unshift!)(B, [ℵ; sort(maxs)])
		println("$ℵ")
	end

	ℵ,ℵ′ = ℵ′, ℵ
	return B, W, ℵ′
end

initIncLoad() = begin ion(); subplot(111,projection="3d"); show(); sleep(.01) end
callbackIncLoad(::Any,::Any,w) = if w!=Void
	plot(w[:,1], w[:,2], w[:,3])
	draw(); sleep(.01)
end

Transient, SteadyState, P	= findCycle(System, .0, 10*rand(3), TransientIterations, TransientStepSize, SteadyStateIterations, SteadyStateStepSize)
y,ω							= prepareCycle(SteadyState, P, fac=1)
Y = fromTraj(y)
V₀ = [Y; ω]
V = matsboNWTN.newton(H, J, V₀, matsboPRED.predCount(48))

B,W,ℵ′ = incrementalLoading(V, .1, 6.0)
B,W,ℵ′ = incrementalLoading(V, .1, 20.0)
# B,W,ℵ′ = incrementalLoading(V, .05, ℵ-1; init=initIncLoad, callback=callbackIncLoad)
# V = matsboNWTN.newton(H, J, V, matsboPRED.predCount(32); init=initPlotResidual, callback=callbackPlotResidual)

# bifurcation diagram
tmpc = map(first, B)
for i in 2:maximum(map(length, B))
	tmp = map(x->i≤length(x)?x[i]:.0, B)
	hold(true)
	plot(tmpc,tmp)
end


# T, SS, ...
begin
	local y = toTraj(V₀[1:end-1])
	ioff(); figure(); ax = subplot(111,projection="3d"); hold(true)
	plot(Transient[:,1], Transient[:,2], Transient[:,3], color="y")
	plot(SteadyState[:,1], SteadyState[:,2], SteadyState[:,3], color="b")
	#plot(cyc[:,1], cyc[:,2], cyc[:,3], color="g")
	plot(y[:,1], y[:,2], y[:,3], color="r")
	ax["scatter"](y[1,1], y[1,2], y[1,3], color="r")
	show(); ion()
end

# Traj before/after cont
figure(); subplot(111,projection="3d"); hold(true)
tmp = toTraj(V[1:end-1]); plot(tmp[:,1], tmp[:,2], tmp[:,3])
tmp = toTraj(W[1:end-1]); plot(tmp[:,1], tmp[:,2], tmp[:,3])

# derivatives
# Jtest = matsboNWTN.broyden(H, v -> matsboNWTN.forwardDifference(H, v))
# Jtest = matsboNWTN.broyden(H, J)
# Jtest = v -> matsboNWTN.forwardDifference(H, v)
# V = matsboNWTN.newton(Hroessler, Jtest, V₀, matsboPRED.predCount(16); init=initPlotResidual, callback=callbackPlotResidual)
