using matsboRK
using matsboPRED
using matsboUTIL

# procedure usually underestimates length of period
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


# cut single cycle of length P*fac from data, resample, shift s.t. X(0)≈0, Fourier transform.
function prepareCycle(data::Array{Float64}, h::Float64, P::Integer; fac::Integer=1)
	local ω, cyc, y

	# cut out singe period
	ω = 2pi / (P*fac * h)
	cyc = data[end-P*fac+2:end,:]

	# resample
	# y = mapslices(V->matsboUTIL.interpolate(V, precision, linspace(1.0, size(cyc,1), samples)), cyc, [1])

	# local I = linspace(1.0, size(cyc,1), Samples)
	# y = mapslices(V->(V[floor(Integer,I)] + V[ceil(Integer,I)])/2, cyc, [1])

	# rotate s.t. x(0)≈0
	# this (↓) is ultimately important! optimization is super sensitive to it...
	y = circshift(cyc,[-findmin(abs(cyc[:,1]))[2]+1])

	return y, ω
end
