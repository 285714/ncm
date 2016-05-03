module matsboUTIL
export mapi, standardize, trim, findPeaks, autocorrelation

mapi(f::Function, v...) = map(f, 1:length(v[1]), v...)
standardize(M::Array) = broadcast((x,y,z)->(x-y)/(z-y), M, minimum(M,[1]), maximum(M,[1]))

function trim(ffront, fback, V)
	local start, stop
	start, stop = 1, size(V,1)
	while ffront(V[start]) start += 1 end
	while ffront(V[stop]) stop -= 1 end
	return start, stop
end


function findPeaks(V)
	local i = 1
	local N = length(V)
	peaks = Integer[]
	while i < length(V)
		# skip non-negative slope area
		while i < N && V[i] <= V[i+1] i += 1 end

		# trace shelf area
		local j = i
		while 1 < j && V[j] == V[j-1] j -= 1 end
		push!(peaks, round(Integer, (j+i)/2))
		i += 1

		# skip non-positive slope area
		while i < length(V) && V[i] >= V[i+1] i += 1 end
	end
	return peaks
end



function autocorrelation(V)
	local tmp = fft(V, [1])
	tmp .*= conj(tmp)
	return real(ifft(tmp, [1]))[1:end÷2,:]
end



end
