# cut single cycle of length P*fac from data, resample, shift s.t. X(0)≈0, Fourier transform.
function prepareCycle(data::Array{Float64}, P::Integer, samples::Integer; fac::Integer=1)
	local ω, cyc, y

	# cut out singe period
	ω = 2pi / (P*fac * SteadyStateStepSize)
	cyc = data[end-P*fac+2:end,:]

	# resample
	y = mapslices(V->matsboUTIL.interpolate(V, 5, linspace(1.0, size(cyc,1), samples)), cyc, [1])

	# rotate s.t. x(0)≈0
	# this (↓) is ultimately important! optimization is super sensitive to it...
	y = circshift(y,[-findmin(abs(y[:,1]))[2]+1])

	return y, ω
end
