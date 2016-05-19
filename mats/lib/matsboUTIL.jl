module matsboUTIL
export mapi, standardize, trim, findPeaks, autocorrelation,
	inputConvert, outputConvert, circconv,
	vectorize, interpolate, interpolateTrigonometric,
	circulantMatrix,
	rows, columns,
	∘, ⊕, ⊗

# map with additional index
# TODO add broadcasting?
mapi(f::Function, v...) = map(f, 1:length(v[1]), v...)

# standardize vector
standardize(V::Array) = broadcast((x,y,z)->(x-y)/(z-y), V, minimum(V,[1]), maximum(V,[1]))

# trim front and/or back of vector V depending on predicates
function trim(ffront, fback, V)
	local start, stop
	start, stop = 1, size(V,1)
	while ffront(V[start]) start += 1 end
	while ffront(V[stop]) stop -= 1 end
	return start:stop
end

# find indices of local maxima, for locally maximal areas return only one representative
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


# autorcorrelate a Vector
# TODO emply rfft?
function autocorrelation(V)
	local tmp = fft(V, [1])
	tmp .*= conj(tmp)
	return real(ifft(tmp, [1]))[1:end÷2,:]
end



# split vector at indices defined in idx (sorted!), apply appropriate contructor
function inputConvert(V, idx, constructor)
    #idx = sort(idx)
    return map((f,i,j)->f(V[i:i+j]), constructor, idx[1:end-1], diff(idx)-1)
end

# deconstruct inputs to vectors and concatenate
function outputConvert(deconstructor, V)
    return reduce(vcat, map((f,v) -> f(v), deconstructor, V))
end
outputConvert(deconstructor, V...) = outputConvert(deconstructor, V)

rows(V) = [ V[i,:] for i in 1:size(V,1) ]
columns(V) = [ V[:,i] for i in 1:size(V,2) ]

# function composition...
# TODO more arguments
∘(f...) = x -> foldr((a,A) -> a(A), x, f)

# parallel composition
⊕(f...) = x -> map(a->f(x), f)


# circular convolution using fft, multiplication in time domain
# TODO employ rfft?
circconv(x,y) = fft(ifft(x) .* ifft(y))
⊗ = circconv

# vectorize a function
vectorize(f) = x -> map(f,x)

# simple periodic (!) Lanczos interpolation
# TODO fix weird modulo thing
function interpolate(V,a::Integer)
	return vectorize() do y
		local sum = zero(V[1])
		local N = length(V)
		for i in floor(y)+(-a+1:a)
			local tmp = y-i
			sum += sinc(tmp)*sinc(tmp/a)*V[((i-1)%N+N)%N + 1]
		end
		return sum
	end
end
interpolate(V, a::Integer, x) = interpolate(V,a)(x)


function circulantMatrix(U)
	local V = reverse(U)
	return reduce(vcat, [ circshift(V,[i])' for i in 1:length(V) ])
end


# returns trigonometric polynomial.
# use with 2a,-2b and divide by 2m+1 to use with rfft coefficients.
function interpolateTrigonometric(a₀, a, b)
	local m = length(a)
	return vectorize() do x
		a₀ + sum(a.*cos(x*(1:m)) + b.*sin(x*(1:m)))
	end
end


end
