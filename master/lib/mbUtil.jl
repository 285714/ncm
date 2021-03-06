__precompile__()

module mbUtil
export mapi, standardize, trim, findPeaks, autocorrelation,
	inputConvert, outputConvert, circconv,
	vectorize,
	circulantMatrix,
	bisection,
	rows, columns,
	∘, ⊕, ⊗,
	@cached

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

function circulantMatrix(C)
	L = length(C)
	return [ C[mod(j-i, L)+1] for i in 1:L, j in 1:L ]
end

function bisection(f,a,b; ϵ=1e-10)
	f(b) < f(a) && ((a,b) = (b,a))
	while true
		c = (a+b)/2
		f(c) < 0 ? a=c : b=c
		(b-a) < ϵ && return (a+b)/2
	end
end


#TODO expose?
varToIter(v...) = return v

function busywait(p)
	t = now()
	l = Base.Dates.Second(p)
	while now()-t < l
		nothing
	end
end

macro cached(m, f)
	#_core = eval(f)
	_name = f.args[1].args[1]
	f.args[1].args[1] = :tmp
	return quote
	$(esc(_name)) = let _cache = Dict(), _core = $f
		function $(esc(_name))(x...)
			#NOTE typechecking is done on caching
			return if haskey(_cache, x)
				_cache[x]
			else
				while length(_cache) ≥ $m
					delete!(_cache, keys(_cache)[1])
				end
				_cache[x] = _core(x...)
			end
		end
	end
	end
end


end #module
