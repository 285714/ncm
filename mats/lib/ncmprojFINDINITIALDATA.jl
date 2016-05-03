# procedure usually underestimates length of period
function findRepresentativeCycle(
	H::Function, t₀::Float64, y₀::Vector,
	TransientIterations::Integer, TransientStepSize::Float64,
	SteadyStateIterations::Integer, SteadyStateStepSize::Float64
	)

	# trace ODE, skipping transient phase
	local data
	local callbackLast(t,y,ft) = data = (t,y)
	local init() = data = Array{Float64}(0,length(y₀))
	local callback(t,y,ft) = data = [data; y']

	if TransientIterations <= 0 data = (t₀,y₀)
	else matsboRK.rk4(H, t₀, y₀, TransientStepSize, matsboPRED.predCount(TransientIterations); callback=callbackLast) end
	matsboRK.rk4(H, data[1], data[2], SteadyStateStepSize, matsboPRED.predCount(SteadyStateIterations); init=init, callback=callback)

	# autocorrelate trajectory, standardize, take minimum
	# 	(and-connect them, only interested in periodicity of ALL components)
	local ac, peaks, P, cyc

	ac = matsboUTIL.autocorrelation(data)
	ac = mapslices(minimum, matsboUTIL.standardize(ac), [2])

	# find peaks in AC, take second highest (not interested in 0-shift)
	peaks = matsboUTIL.findPeaks(ac)[2:end-1]
	P = peaks[findmax(ac[peaks])[2]]
	cyc = data[end-P+1:end,:]

	return data, cyc
end
