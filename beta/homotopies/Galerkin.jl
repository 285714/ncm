#NOTE Galerkin defined this way allows only for numerical derivation.

type Galerkin <: Homotopy
	f::Function
	d::Int
end

#NOTE untested

function H(G::Galerkin)
	local m = ((length(V)-2)÷G.d - 1) ÷ 2
	return function H(V::Vector{Float64})
		local C = reduce(hcat, [ complex(V[(i-1)*(2m+1)+(1:m+1)], [0, V[(i-1)*(2m+1)+(m+1)+(1:m+1)]]) for i in 1:G.d ])
		C = rfft(mapslices(G.f, irfft(C, 2m, [1]), [2]), [1])
		return reduce(vcat, [ [real(C[i]); imag(C[i][2:end])] for i in 1:G.d ])
	end
end

function J(G::Galerkin)
	error("not defined")
end

function V(G::Galerkin)
	error("not defined")
end

function C(G::Galerkin)
	error("not defined")
end
