# stable periodic trajectories for ρ ∈ {350, 160, 100.5, 9.65}

global σ=10, β=8/3
lorenz(t,v) = [
	σ*(v[2]-v[1])
	v[1]*(v[4]-v[3])-v[2]
	v[1]*v[2]-β*v[3]
	]
lorenz(v) = lorenz(0, v)

function fToH(f)
	return function H(V::Vector{Float64})
		m = length(V-5)÷6
		anchor = sum(V[1:m+1])
		ω,ρ = V[end-1:end]
		V = reshape(V[1:end-2],2m+1,3)
		V = complex(V[1:m+1,:], [zeros(3)'; V[m+2:end,:]])

		defect = rfft(mapslices(x->f([x;ρ]), irfft(V, 2m, [1]), [2]), [1]) - im*ω*(0:m) .* V
		defect = vec([ real(defect); imag(defect)[2:end,:] ])

		return [defect; anchor]
	end
end

##

f = lorenz
H = fToH(lorenz)
J(V) = forwardDifference(H, V)
