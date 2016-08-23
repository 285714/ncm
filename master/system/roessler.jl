const α = .1
const β = .1
roessler(t,v) = [
	-v[2] - v[3]
	v[1] + α*v[2]
	β + v[3]*(v[1]-v[4])
	]

roessler′(t,v) = [
	0		-1		-1				0
	1		α		0				0
	v[3]	0		(v[1]-v[4])	 	-v[3]
	]

function Hroessler(V::Vector{Float64})
	local m = length(V-5)÷6
	local X₀,Xᵣ,Xᵢ,Y₀,Yᵣ,Yᵢ,Z₀,Zᵣ,Zᵢ,ω,ℵ

	ω,ℵ = V[end-1:end]
	V = reshape(V[1:end-2],2m+1,3)
	X₀,Y₀,Z₀ = V[1,1], V[1,2], V[1,3]
	Xᵣ,Yᵣ,Zᵣ = V[2:2+m-1,1], V[2:2+m-1,2], V[2:2+m-1,3]
	Xᵢ,Yᵢ,Zᵢ = V[2+m:end,1], V[2+m:end,2], V[2+m:end,3]

	local S₀,Sᵣ,Sᵢ
	local tmp = rfft( irfft(complex([X₀;Xᵣ], [.0;Xᵢ]), 2m+1) .* irfft(complex([Z₀;Zᵣ], [.0;Zᵢ]), 2m+1) )
	S₀,Sᵣ,Sᵢ = real(tmp[1]), real(tmp[2:end]), imag(tmp[2:end])

	local D = ω*(1:m)

	rtn = [
		#A
		-Y₀-Z₀
		-Yᵣ-Zᵣ+D.*Xᵢ
		-Yᵢ-Zᵢ-D.*Xᵣ
		#B
		X₀+α*Y₀
		Xᵣ+α*Yᵣ+D.*Yᵢ
		Xᵢ+α*Yᵢ-D.*Yᵣ
		#C
		-ℵ*Z₀+S₀+(2m+1)*β
		-ℵ*Zᵣ+Sᵣ+D.*Zᵢ
		-ℵ*Zᵢ+Sᵢ-D.*Zᵣ

		X₀+sum(Xᵣ)
	]

	return rtn
end



function Jroessler(V)
	local m = length(V-5)÷6
	local X₀,Xᵣ,Xᵢ,Y₀,Yᵣ,Yᵢ,Z₀,Zᵣ,Zᵢ,ω,ℵ

	ω,ℵ = V[end-1:end]
	V = reshape(V[1:end-2],2m+1,3)
	X₀,Y₀,Z₀ = V[1,1], V[1,2], V[1,3]
	Xᵣ,Yᵣ,Zᵣ = V[2:2+m-1,1], V[2:2+m-1,2], V[2:2+m-1,3]
	Xᵢ,Yᵢ,Zᵢ = V[2+m:end,1], V[2+m:end,2], V[2+m:end,3]

	local D = ω*(1:m)

	# Jacobian of circular convolution of coefficients of real functions in appropriate format.
	# ∂/∂X cconv(X,Y) = ∂/∂X F(x⋅y) = rCCD(Y)
	rCCD(C) = rCCD(real(C[1]), real(C[2:end]), imag(C[2:end]))
	function rCCD(V₀,Vᵣ,Vᵢ)
		local I1,I2,Wᵣ,Wᵢ
		I1 = [ mod(i-j, 2m+1)+1 for i in 0:m, j in 0:m ]
		I2 = [ mod(i+j, 2m+1)+1 for i in 0:m, j in 0:m ]
		Wᵣ = [V₀; Vᵣ; Vᵣ[end:-1:1]]
		Wᵢ = [.0; Vᵢ; -Vᵢ[end:-1:1]]

		local rtn =  [
			(Wᵣ[I1]+Wᵣ[I2])					(-Wᵢ[I1]+Wᵢ[I2])[:,2:end]
			(Wᵢ[I1]+Wᵢ[I2])[2:end,:]		(Wᵣ[I1]-Wᵣ[I2])[2:end,2:end]
		] / (2m+1)
		rtn[:,1] /= 2.0
		return rtn
	end

	local AbyX, AbyY, AbyZ, Abyω, Abyℵ,
		BbyX, BbyY, BbyZ, Bbyω, Bbyℵ,
		CbyX, CbyY, CbyZ, Cbyω, Cbyℵ,
		LbyX, LbyY, LbyZ, Lbyω, Lbyℵ

	AbyX = [
		.0					zeros(1,m)			zeros(1,m)
		zeros(m)			zeros(m,m)			diagm(D)
		zeros(m)			diagm(-D)			zeros(m,m)
	]

	AbyY = AbyZ = -eye(2m+1)

	BbyX = eye(2m+1)

	BbyY = [
		α					zeros(1,m)			zeros(1,m)
		zeros(m)			α*eye(m)			diagm(D)
		zeros(m)			diagm(-D)			α*eye(m)
	]

	BbyZ = zeros(2m+1,2m+1)

	CbyX = rCCD(Z₀,Zᵣ,Zᵢ)

	CbyY = zeros(2m+1, 2m+1)

	CbyZ = rCCD(X₀,Xᵣ,Xᵢ) + [
		-ℵ					zeros(1,m)			zeros(1,m)
		zeros(m)			-ℵ*eye(m,m)			diagm(D)
		zeros(m)			diagm(-D)			-ℵ*eye(m,m)
	]

	LbyX = [ 1.0			ones(1,m)			zeros(1,m) ]
	LbyY = zeros(1,2m+1)
	LbyZ = zeros(1,2m+1)

	Abyω = [
		.0
		(1:m).*Xᵢ
		-(1:m).*Xᵣ
	]

	Bbyω = [
		.0
		(1:m).*Yᵢ
		-(1:m).*Yᵣ
	]

	Cbyω = [
		.0
		(1:m).*Zᵢ
		-(1:m).*Zᵣ
	]

	Lbyω = .0

	Abyℵ = Bbyℵ = zeros(2m+1)
	Cbyℵ = [-Z₀; -Zᵣ; -Zᵢ]
	Lbyℵ = .0

	return [
		AbyX AbyY AbyZ Abyω Abyℵ
		BbyX BbyY BbyZ Bbyω Bbyℵ
		CbyX CbyY CbyZ Cbyω Cbyℵ
		LbyX LbyY LbyZ Lbyω Lbyℵ
	]
end



function roesslerProjectionInternal(V)
	W,ω,ρ = unwrap(V)
	m = size(W,1)-1
	ftmp,N = interps(W), 2m

	const p = 0.08872374069251765
	T = linspace(-pi+p,pi+p,N)

	pA(t) = ftmp[1](t)
	# pB(t) = ftmp[1](t) + ftmp[2](t) - 2d
	t₀ = T[end]
	a₀ = pA(t₀)
	# b₀ = pB(t₀)

	ts,fs,ns = [],[],[]
	for t in T
		tmp = pA(t)
		if a₀ ≤ 0 < tmp
			push!(ts, bisection(pA, t₀, t))
			push!(fs, map(g->g(ts[end]), ftmp))
			push!(ns, norm(fs[end]))
		end
		a₀ = tmp

		# tmp = pB(t)
		# if b₀ ≤ 0 < tmp
		# 	push!(ts, bisection(pB, t₀, t))
		# 	push!(fs, map(g->g(ts[end]), ftmp))
		# 	push!(ns, -norm(fs[end]-B))
		# end
		# b₀ = tmp

		t₀ = t
	end

	global roesslerfs = fs, roesslerstart = map(f->f(.0), ftmp)

	return ts, fs, sort(ns)
end


roesslerProjection(v...) = roesslerProjectionInternal(v...)[3]
function cbPlotSolution(v)
	tmp = reduce(hcat, roesslerfs)'
	scatter3D(tmp[:,1],tmp[:,2],tmp[:,3])
	scatter3D(roesslerstart[1], roesslerstart[2], roesslerstart[3], c="r")
end
