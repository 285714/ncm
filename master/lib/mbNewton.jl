__precompile__()

module mbNewton
export newton, centralDifference, forwardDifference, broyden, bisection

# base function
"""
    newton(homotopy, jacobian, v₀, pred[, init, callback, useOpt])
Newton's method. For example: `v1 = newton(H, J, v0, predEps(0.001))`.
"""
function newton(H::Function, J::Function, v₀::Vector{Float64}, pred::Function;
	init=Void, callback=Void, useOpt::Bool=false)

	local v = deepcopy(v₀)
	local tmpH, tmpJ
	tmpH, tmpJ = H(v), J(v)

	init≠Void && init()
	callback≠Void && callback(v, tmpH, tmpJ)

	if useOpt
		local optv, optH, tmpH2
		optv, optH = v, norm(tmpH)
	end

	while pred(tmpH)
		# v -= size(tmpJ,1)==size(tmpJ,2) ? tmpJ\tmpH : tmpJ'*inv(tmpJ*tmpJ')*tmpH
		v -= tmpJ\tmpH

		tmpH, tmpJ = H(v), J(v)

		callback≠Void && callback(v, tmpH, tmpJ)

		if useOpt && (tmpH2 = norm(tmpH)) < optH
			optv, optH = v, tmpH2
		end
	end

  return useOpt ? optv : v
end


"""
    centralDifference(homotopy, v, epsilon)
	centralDifference(homotopy, epsilon)
Same as `forwardDifference`, but with the symmetric difference quotient.
"""
function centralDifference(H::Function, v::Vector{Float64}; ϵ::Float64=1e-4)
	local J = cell(length(v))
	local w = deepcopy(v)
	for i in 1:length(v)
		w[i] = v[i]+ϵ; 	J[i] = H(w)
		w[i] = v[i]-ϵ; 	J[i] -= H(w)
	 	J[i] /= 2ϵ;		w[i] = v[i]
	end
	return reduce(hcat, J)
end
centralDifference(H::Function; ϵ=1e-4) = v -> centralDifference(H, v; ϵ)

"""
    forwardDifference(homotopy, v, epsilon)
	forwardDifference(homotopy, epsilon)
Returns the difference quotient in `v` or an approximate jacobian using this method.
"""
function forwardDifference(H::Function, v::Vector{Float64}; ϵ::Float64=1e-4)
	local J = cell(length(v))
	local w = deepcopy(v)
	local tmpH = H(v)
	for i in 1:length(v)
	  w[i] = v[i]+ϵ;	J[i] = (H(w) - tmpH) ./ ϵ
	  w[i] = v[i]
	end
	return reduce(hcat, J)
end
forwardDifference(H::Function; ϵ=1e-4) = v -> forwardDifference(H, v; ϵ=ϵ)


# WARNING inefficient, H is evaluated which is not necessary (at least in Newton context...
#		asymptotic behavior is unchanged though) however, this allows to construct a pseudo-
#		Jacobian with the usual one-parameter-signature, instead of requiring to adapt any functions.
# can be countered through rudimentary caching in H
# WARNING mind that this function has an internal state depending on the preceding
#		evaluations. the 'Jacobian' at a point is thus not well defined.
# TODO sequence of evaluation not optimal
"""
    broyden(homotopy, jacobian)
Creates an approximate jacobian based on broyden-updates and maintainance of an internal
state. Consider this for performance improvements.
"""
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

		norm(dv)>1e-10 && (J₀ += (dH-J₀*dv) ./ norm(dv)^2 * dv')

		return J₀
	end
end




end
