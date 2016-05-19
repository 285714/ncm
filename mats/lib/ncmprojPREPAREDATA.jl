# cut single cycle of length P*fac from data, resample, shift s.t. X(0)≈0, Fourier transform.
function prepareCycle(data::Array{Float64}, P::Integer; fac::Integer=1)
	local ω, cyc, y

	# cut out singe period
	ω = 2pi / (P*fac * SteadyStateStepSize)
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
