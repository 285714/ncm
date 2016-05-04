module matsboNWTN
export newton, centralDifference, forwardDifference

# base function
function newton(H::Function, J::Function, v₀::Vector{Float64}, pred::Function; init::Function=()->Void, callback::Function=(v...)->Void)
    local v = deepcopy(v₀)
    local tmpH, tmpJ
    tmpH, tmpJ = H(v), J(v)

    init()
    callback(v, tmpH, tmpJ)

    while pred(tmpH)
        v -= tmpJ \ tmpH
        tmpH, tmpJ = H(v), J(v)
        callback(v, tmpH, tmpJ)
    end

    return v
end


# numerical differentiation methods
function centralDifference(H::Function, v::Vector{Float64}; ϵ::Float64=.001, fac::Float64=1.0)
    local J = cell(length(v))
		local w = deepcopy(v)
		local tmp = fac/2ϵ
    for i in 1:length(v)
				w[i] = v[i]+ϵ; 	J[i] = H(w)
				w[i] = v[i]-ϵ; 	J[i] -= H(w)
        J[i] *= tmp; 		w[i] = v[i]
    end
    return reduce(hcat, J)
end

function forwardDifference(H::Function, v::Vector{Float64}; ϵ::Float64=.001, fac::Float64=1.0)
    local J = cell(length(v))
    local w = deepcopy(v)
    local tmpH = H(v)
		local tmp = fac/ϵ
    for i in 1:length(v)
        w[i] = v[i]+ϵ;	J[i] = (H(w) - tmpH) * tmp
        w[i] = v[i]
    end
    return reduce(hcat, J)
end

end
