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

function f′ToJ(f′)
	return function J(V::Vector{Float64})
		m = length(V-5)÷6
		ω,ρ = V[end-1:end]
		V = reshape(V[1:end-2],2m+1,3)
		V = complex(V[1:m+1,:], [zeros(3)'; V[m+2:end,:]])

		T = irfft(V, 2m, [1])
		M = Array{Float64}(size(T,1), size(T,2), size(T,2)+1) #n, comp, deriv
		for i in 1:size(T,1) M[i,:,:] = f′([T[i,:]'; ρ]) end

		M′ = [ rCCD(rfft(M[:,i,j])) for i in 1:size(T,2), j in 1:size(T,2) ]

		for i in 1:size(T,2)
			D = diagm(ω*(0:m))
			M′[i,i][m+2:end,1:m+1] -= D[2:end,:]
			M′[i,i][1:m+1, m+2:end] += D[:,2:end]
		end

		M′ = reducedim(vcat, M′, [1], Array{Float64}(0,size(M′[1,1],2)))
		M′ = reducedim(hcat, M′, [2], Array{Float64}(size(M′[1,1],1),0))[1]

		ddω = -im*(0:m).*V
		ddω = vec([real(ddω); imag(ddω)[2:end,:]])

		ddρ = reduce(hcat, [ rfft(M[:,i,size(T,2)+1]) for i in 1:size(T,2) ])
		ddρ = vec([real(ddρ); imag(ddρ)[2:end,:]])

		return [
			M′ ddω ddρ
			[ones(1,m+1) zeros(1,5m+2)] 0 0
		]
	end
end




# Jacobian of circular convolution of coefficients of real functions in appropriate format.
# ∂/∂X cconv(X,Y) = ∂/∂X F(x⋅y) = rCCD(Y)
rCCD(C) = rCCD(real(C[1]), real(C[2:end]), imag(C[2:end]))
function rCCD(V₀,Vᵣ,Vᵢ)
	m = length(Vᵣ)
	I1 = [ mod(i-j, 2m+1)+1 for i in 0:m, j in 0:m ]
	I2 = [ mod(i+j, 2m+1)+1 for i in 0:m, j in 0:m ]
	Wᵣ = [V₀; Vᵣ; Vᵣ[end:-1:1]]
	Wᵢ = [.0; Vᵢ; -Vᵢ[end:-1:1]]

	local rtn =  [
		(Wᵣ[I1]+Wᵣ[I2])					(-Wᵢ[I1]+Wᵢ[I2])[:,2:end]
		(Wᵢ[I1]+Wᵢ[I2])[2:end,:]		(Wᵣ[I1]-Wᵣ[I2])[2:end,2:end]
	] / (2m)
	rtn[:,1] /= 2
	return rtn
end
