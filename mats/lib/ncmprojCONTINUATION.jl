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
		W = matsboNWTN.newton(H, J, W, predCount(8)∧predEps(1e-4))
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




function PC(H, J, V₀, pred; init=Void, callback=Void, h₀=1.0, κ₀=.5, δ₀=1.0, α₀=2pi*10/360, dir=true)
	local V,h; V,h = deepcopy(V₀),h₀

	init≠Void && init()
	callback≠Void && callback(V)

	local T′ = Void
	while pred(V)
		# P-step
		local T = matsboUTIL.tangent(J(V))
		local W = V + (dir?h:-h) * T

		# step-length adaption
		local κ,δ,α

		T′==Void && (T′=T)
		α = acos(clamp(dot(T′, T), -1.0, 1.0))
		T′=T

		local X = Any[]
		W = matsboNWTN.newton(H, J, W, predCount(2); callback=(W,H,J)->push!(X,W))
		δ = norm(X[1]-X[2])
		κ = norm(X[3]-X[2]) / δ

		local f = clamp(max(sqrt(κ/κ₀), sqrt(δ/δ₀), α/α₀), .5, 2.0)
		h /= f

		println("$((sqrt(κ/κ₀), sqrt(δ/δ₀), α/α₀))")

		# C-step
		f==2.0 && continue
		V = matsboNWTN.newton(H, J, W, predEps(1e-2))

		# println("$h")

		callback≠Void && callback(V)
	end

	return V
end
