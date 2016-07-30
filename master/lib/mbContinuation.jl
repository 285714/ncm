__precompile__()

module mbContinuation
export PC, tangent, incrementalLoading

using mbPred, mbNewton

function tangent(A)
	local t = [A;ones(1,size(A,2))] \ [zeros(size(A,1)); 1.0]
	return flipsign(t / norm(t), det([A;t']))
end


function PC(H, J, V₀, pred; init=Void, callback=Void, predC=predEps(1e-10), h₀=1.0, κ₀=.5, δ₀=1.0, α₀=pi/36, dir=true)
	local V,h; V,h = deepcopy(V₀),h₀
	V = mbNewton.newton(H, J, V, predC)

	init≠Void && init()
	callback≠Void && callback(V)

	local T′ = Void
	while pred(V)
		# P-step
		local T = tangent(J(V))
		local W = V + (dir?h:-h) * T

		# step-length adaption
		local κ,δ,α

		T′==Void && (T′=T)
		α = acos(clamp(dot(T′, T), -1.0, 1.0))
		T′=T

		α>.95pi && break #TODO fix... bifurcation point detection...

		local X = Any[]
		W = mbNewton.newton(H, J, W, predCount(2); callback=(W,H,J)->push!(X,W))
		δ = norm(X[1]-X[2])
		κ = norm(X[3]-X[2]) / δ

		local f = clamp(max(sqrt(κ/κ₀), sqrt(δ/δ₀), α/α₀), .5, 2.0)
		h /= f

		# C-step
		f≥2.0 && continue
		V = mbNewton.newton(H, J, W, predC)

		callback≠Void && callback(V)
	end

	return V
end


#NOTE untested!
function incrementalLoading(H, J, V₀, R; init=Void, callback=Void, predC=predEps(1e-10))
	local V = deepcopy(V₀)

	init≠Void && init()
	callback≠Void && callback(V)

	for ℵ in R
		V = mbNewton.newton(H, J, V, predC)
		callback≠Void && callback(V)
	end

	return V
end



end #module
