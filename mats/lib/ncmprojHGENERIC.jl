# WARNING wip!

function Hgeneric(f::Function, dims::Integer)
	function rToC(V)
		local N = (length(V)-1)÷dims
		local c(x) = [ x[1]+0im; complex(x[2:end÷2+1], x[end÷2+2:end]) ]
		local tmp = matsboUTIL.inputConvert(V, [1+N*(0:dims); dims*N+2], [collect(repeated(c,dims)); identity])
		return reduce(hcat, tmp[1:end-1]), tmp[end]
	end

	function cToR(V)
		local d(x) = [real(x[1]); real(x[2:end]); imag(x[2:end])]
		return matsboUTIL.outputConvert([collect(repeated(d,dims)); identity], V)
	end

	return function H(V)
		local N,C,ω
		N = (length(V)-1)÷dims
		C, ω = rToC(V)

		local tmp = rfft(mapslices(f, irfft(C, N, [1]), [2]), [1])
		tmp = broadcast((a,b,c)->a-b*c, tmp, im*ω.*(0:N÷2), C)

		return cToR([ matsboUTIL.columns(tmp); reduce((X,x)->X+real(x), .0, C[:,1]) ])
	end
end
