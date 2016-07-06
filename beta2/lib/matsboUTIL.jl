__precompile__()

module matsboUTIL
export mapi, standardize, trim, findPeaks, autocorrelation,
	inputConvert, outputConvert, circconv,
	vectorize,
	circulantMatrix,
	bisection,
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

# TODO inefficient
function circulantMatrix(U)
	local V = reverse(U)
	return reduce(vcat, [ circshift(V,[i])' for i in 1:length(V) ])
end

function bisection(f,a,b; ϵ=1e-10)
	b < a && ((a,b) = (b,a))
	while norm(a-b) > ϵ
		local c = (a+b)/2
		f(c) < 0 ? a=c : b=c
	end
	return a
end

# prefer the compact solution, plus stays in sync.
# function bisection(f, x₋, x₊; ϵ=1e-10)
# 	while true
# 		x₀ = (x₋ + x₊) / 2
# 		if f(x₀) == 0 || abs(x₋ - x₊) < ϵ
# 			return x₀
# 		elseif sign(f(x₀)) * sign(f(x₋)) < 0
# 			x₊ = x₀
# 		else
# 			x₋ = x₀
# 		end
# 	end
# end


end
