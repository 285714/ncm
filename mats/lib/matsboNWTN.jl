module matsboNWTN
export newton, centralDifference, forwardDifference, broyden

# base function
function newton(H::Function, J::Function, v₀::Vector{Float64}, pred::Function;
	init::Function=()->Void, callback::Function=(v...)->Void, useOpt::Bool=true)

	local v = deepcopy(v₀)
	local tmpH, tmpJ
	tmpH, tmpJ = H(v), J(v)

	init()
	callback(v, tmpH, tmpJ)

	if useOpt
		local optv, optH, tmpH2
		optv, optH = v, norm(tmpH)
	end

	while pred(tmpH)
	  v -= tmpJ \ tmpH
	  tmpH, tmpJ = H(v), J(v)

	  callback(v, tmpH, tmpJ)

		if useOpt && (tmpH2 = norm(tmpH)) < optH
			optv, optH = v, tmpH2
		end
	end

  return useOpt ? optv : v
end


# numerical differentiation methods
function centralDifference(H::Function, v::Vector{Float64}; ϵ::Float64=1e-4)
	local J = cell(length(v))
	local w = deepcopy(v)
	for i in 1:length(v)
		w[i] = v[i]+ϵ; 	J[i] = H(w)
		w[i] = v[i]-ϵ; 	J[i] -= H(w)
	  J[i] /= 2ϵ; 		w[i] = v[i]
	end
	return reduce(hcat, J)
end

function forwardDifference(H::Function, v::Vector{Float64}; ϵ::Float64=1e-4)
	local J = cell(length(v))
	local w = deepcopy(v)
	local tmpH = H(v)
	for i in 1:length(v)
	  w[i] = v[i]+ϵ;	J[i] = (H(w) - tmpH) / ϵ
	  w[i] = v[i]
	end
	return reduce(hcat, J)
end


# WARNING inefficient, H is evaluated which is not necessary (at least in Newton context...
#		asymptotic behavior is unchanged though) however, this allows to construct a pseudo-
#		Jacobian with the usual one-parameter-signature, instead of requiring to adapt any functions.
# WARNING mind that this function has an internal state depending on the preceding
#		evaluations. the 'Jacobian' at a point is thus not well defined.
# TODO sequence of evaluation not optimal
function broyden(H, J)
	local H₀,J₀,v₀							# closure
	local init = true

	return function (v₁)
		if init
			init = false
			H₀,J₀,v₀ = H(v₁),J(v₁),v₁
			return J₀
		end

		local H₁, dv, dH
		H₁ = H(v₁) 								# ← hmm...
		dv, dH = v₁-v₀, H₁-H₀
		v₀, H₀ = v₁, H₁

		return J₀ += (dH-J₀*dv) / norm(dv)^2 * dv'
	end
end




end
