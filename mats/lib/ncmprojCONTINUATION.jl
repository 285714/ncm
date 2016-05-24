function toTrajInterp(V)
	global ExternSamples
	local m = (length(V)-Dimensions)÷(2*Dimensions)
	local tmp = reshape(V, 2m+1, Dimensions)
	local P = map(1:Dimensions) do i
		x -> interpolateTrigonometric(tmp[1,i], 2tmp[2:2+m-1,i], -2tmp[2+m:end,i])(x) / (2m+1)
	end
	return P
end


function toTraj(V)
	global ExternSamples
	local P = toTrajInterp(V)
	local x = linspace(.0, 2pi, ExternSamples)
	return reduce(hcat, map(f->f(x), P))
end

function fromTraj(v)
	global InterpPrec, InternSamples
	local y = mapslices(v->interpolateLanczos(v, InterpPrec, linspace(1.0, size(v,1), InternSamples)), v, [1])
	local Y = rfft(y, [1])
	return vec(mapslices(y -> [real(y); imag(y[2:end])], Y, [1]))
end
