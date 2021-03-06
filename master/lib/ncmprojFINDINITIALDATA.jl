using mbRK
using mbPred
using mbUtil

# procedure usually underestimates length of period
"""
    findCycle(H, t0, y0, transientIterations, transientStepSize,
	    steadyStateIterations, steadyStateStepSize)
Returns all points from the `transientIterations`-times numerical integration of `H` with
fixed step-size `transientStepSize`, and all points from the `steadyStateIterations`-times
integration with fixed step-size `steadyStateStepSize` (assuming it hit the steady state
part), and the periodicity of the found cycle.
"""
function findCycle(
	H::Function, t₀::Float64, y₀::Vector,
	TransientIterations::Integer, TransientStepSize::Float64,
	SteadyStateIterations::Integer, SteadyStateStepSize::Float64
	)

	local dataTransient, dataSteadyState
	dataTransient, dataSteadyState = Array{Float64}(0,length(y₀)), Array{Float64}(0,length(y₀))
	function callbackTransient(t,y,ft) dataTransient = [dataTransient; y'] end
	function callbackSteadyState(t,y,ft) dataSteadyState = [dataSteadyState; y'] end

	# transient part
	if TransientIterations <= 0
		dataTransient = y₀'
	else
		# rk1(H, t₀, y₀, TransientStepSize/10, predCount(TransientIterations); callback=callbackTransient)
		rk4(H, t₀, y₀, TransientStepSize, predCount(TransientIterations); callback=callbackTransient)
	end

	# (hopefully) steady state part
	local t₁, y₁
	t₁, y₁ = TransientIterations*TransientStepSize, squeeze(dataTransient[end,:]', (2))
	# rk1(H, t₁, y₁, SteadyStateStepSize/10, predCount(SteadyStateIterations); callback=callbackSteadyState)
	rk4(H, t₁, y₁, SteadyStateStepSize, predCount(SteadyStateIterations); callback=callbackSteadyState)

	# autocorrelate trajectory, standardize, take minimum
	# 	(and-connect them, only interested in periodicity of ALL components)
	local ac, peaks, P, cyc

	ac = autocorrelation(dataSteadyState)
	ac = mapslices(minimum, standardize(ac), [2])

	# find peaks in AC, take second highest (not interested in 0-shift)
	peaks = findPeaks(ac)[2:end-1]
	P = peaks[findmax(ac[peaks])[2]]

	return dataTransient, dataSteadyState, P
end


"""
    prepareCycle(data, h, P[, fac])
cut single cycle of length `P*fac` from data, resample, shift s.t. `X(0)≈0`, Fourier transform.
"""
function prepareCycle(data::Array{Float64}, h::Float64, P::Integer; fac::Integer=1)
	local ω, cyc, y

	# cut out singe period
	ω = 2pi / (P*fac * h)
	cyc = data[end-P*fac+2:end,:]

	# resample
	# y = mapslices(V->mbUtil.interpolate(V, precision, linspace(1.0, size(cyc,1), samples)), cyc, [1])

	# local I = linspace(1.0, size(cyc,1), Samples)
	# y = mapslices(V->(V[floor(Integer,I)] + V[ceil(Integer,I)])/2, cyc, [1])

	# rotate s.t. x(0)≈0
	# this (↓) is ultimately important! optimization is super sensitive to it...
	y = circshift(cyc,[-findmin(abs(cyc[:,1]))[2]+1])

	return y, ω
end


"""
    findCyclePoincare(F, y[, plane, clusterRating, nIntersections,
	    maxCycles, sampleSize, transientIterations, transientStepSize,
		steadyStateStepSize])
extracts a single cycle of the steady state of ode `F` using poincare cuts through
the `plane`.
"""
function findCyclePoincare(
		F :: Function,
		y₀ :: Vector;
		plane :: Function = first,	# x ∈ plane ⇔ plane(x) = 0
		clusterRating :: Function = x -> var(x) / length(x),	# use penalty for more points, i.e., normalize
		nIntersections :: Int = 120,
		maxCycles :: Int = 30,
		sampleSize :: Int = 512,
		transientIterations :: Int = 2000,
		transientStepSize :: Float64 = 0.1,
		steadyStateStepSize :: Float64 = 0.1)

	dim = length(y₀)
	flagIntersected = false
	y⁻ = y₀
	y₁ = ode(F, y₀, transientStepSize, function (y, i, t, V)
		plane(y⁻) ≤ 0 ≤ plane(y) && (flagIntersected = true)
		y⁻ = y
		return i < transientIterations
	end)

	flagIntersected || error("No intersection in transient phase.")

	# generate intersections with plane
	h = steadyStateStepSize
	y⁻ = y₁
	intersections = zeros(dim + 1, 0)
	ode(F, y₁, h, (y, i, t, V) -> begin
		if (sign(plane(y)) > 0 &&  sign(plane(y⁻)) <= 0)
			# find h′ ∈ [0,h] s.t. plane(y⁻ + h′⋅V(h′,y⁻)) = 0
			step(h′) = y⁻ + h′*V(h′,y⁻)
			h′ = mbUtil.bisection(h′ -> plane(step(h′)), 0, h, ϵ=1e-4)
			y′ = step(h′)
			intersections = [intersections [y′; t+h′]]
		end
		y⁻ = y
		i > 100*transientStepSize/steadyStateStepSize*transientIterations && error("No intersection in steady state phase.")
		nIntersections > size(intersections, 2)
	end)

	#scatter(intersections[2,:], intersections[3,:])

	# find number of cycles s.t. rating is minimum
	n = size(intersections, 2)
	ratings = map(1:maxCycles) do i
		mean( [ sum(var(intersections[1:end-1, j:i:n], [2])) for j in 1:i ] )
	end

	mn₁,m₁ = minimum(ratings), mean(ratings)
	period = findfirst(1:maxCycles) do i
		X = ratings[i:i:end]
		# mn₂,m₂ = minimum(ratings), mean(ratings)
		std(X)/mean(X) < .5 && all(X) do x
			abs(x-mn₁)/m₁ < .5
		end
	end

	if period == 0
		period = 1
		warn("A valid period could not be found.")
	end

	# cut out steady state trajectory
	y₂ = intersections[1:dim,1]
	T = intersections[dim+1,1+period] - intersections[dim+1,1]
	dataSteadyState = zeros(0,dim)
	ode(F, y₂, steadyStateStepSize, (y, i, t, V) -> begin
		dataSteadyState = [dataSteadyState; y']
		return t < T
	end)

	dataSteadyState, period, 2π/T
end


# externalize or merge with rk4
function ode(f, y, h, pred)
	RK  = [1/2 0 0 0; 0 1/2 0 0; 0 0 1 0; 1/6 2/6 2/6 1/6]; # classic
	s   = size(RK, 1);
	dim = length(y);
	t   = 0

	function V(h, y)
		k = zeros(dim, s);

		k[:,1] = f(t, y);
		for j=1:s-1
			k[:,j+1] = f( t + h * sum(RK[j,:]), y + h * k * RK[j,:]' );
		end

		k * RK[s,:]';
	end

	i = 0
	while pred(y, i, t, V) !== false
		y = y + h * V(h, y);
		t += h
		i += 1
	end

	y
end


# y = rand(3)
# T = zeros(0,3)
# y = ode((t,v)->f(t,[v;350]), y, .01, (y, i, t, V) -> begin
# 	global T = [T; y']
# 	return i < 1000
# end); println(y)
#
# f = figure()
# gca(projection="3d")
# plot(T[:,1], T[:,2], T[:,3])
# f[:show]()
