function toTraj(V)
	global ExternSamples
	local m = (length(V)-Dimensions)÷(2*Dimensions)
	local tmp = reshape(V, 2m+1, Dimensions)
	local P = map(1:Dimensions) do i
		x -> matsboUTIL.interpolateTrigonometric(tmp[1,i], 2tmp[2:2+m-1,i], -2tmp[2+m:end,i])(x) / (2m+1)
	end
	local x = linspace(.0, 2pi, ExternSamples)
	return reduce(hcat, map(f->f(x), P))
end

function fromTraj(v)
	global InterpPrec, InternSamples
	local y = mapslices(v->matsboUTIL.interpolate(v, InterpPrec, linspace(1.0, size(v,1), InternSamples)), v, [1])
	local Y = rfft(y, [1])
	return vec(mapslices(y -> [real(y); imag(y[2:end])], Y, [1]))
end

function incrementalLoading(V, Δℵ, ℶ; init=Void, callback=Void)
	global ℵ
	local B = Array{Float64}[]
	local W = deepcopy(V)
	local ℵ′ = ℵ

	init≠Void && init()
	callback≠Void && callback(ℵ, W, Void)

	for ℵ′′ in ℵ:(ℵ≤ℶ?Δℵ:-Δℵ):ℶ 	# =)
		ℵ = ℵ′′
		W = matsboNWTN.newton(H, J, W, predCount(128)∧predEps(.001))
		local w = toTraj(W[1:end-1])
		callback≠Void&& callback(ℵ, W, w)
		local maxs = Float64[]
		for i in 1:size(w,1)-1
			w[i,1] ≤ 0 && w[i+1,1] > 0 && push!(maxs, norm(w[i,:]))
		end
		w[end,1] ≤ 0 && w[1,1] > 0 && push!(maxs, norm(w[end,:]))
		(ℵ′<ℶ ? push! : unshift!)(B, [ℵ; sort(maxs)])
		println("$ℵ") #DEBUG
	end

	ℵ,ℵ′ = ℵ′, ℵ
	return B, W, ℵ′
end




function PC(H,J,V₀,pred)
	local h,V
	h = .001
	V = deepcopy(V₀)

	while pred
		V += h*matsboUTIL.tangent(J(V))
		V = matsboNWTN.newton(H,J,V,predCount(128)∧predEps(.001))
	end

	return V
end
