a,b = .1,.1
roessler(v::Vector) = [-v[2] - v[3]; v[1] + a*v[2]; b + v[3]*(v[1]-ℵ)]
roessler(t,v) = roessler(v)

function Hroessler(V::Vector{Float64})
	local m = length(V-4)÷6
	local ω,X₀,Xᵣ,Xᵢ,Y₀,Yᵣ,Yᵢ,Z₀,Zᵣ,Zᵢ

	ω = V[end]
	V = reshape(V[1:end-1],2m+1,3)
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
		X₀+a*Y₀
		Xᵣ+a*Yᵣ+D.*Yᵢ
		Xᵢ+a*Yᵢ-D.*Yᵣ
		#C
		(2m+1)*b-ℵ*Z₀+S₀
		-ℵ*Zᵣ+Sᵣ+D.*Zᵢ
		-ℵ*Zᵢ+Sᵢ-D.*Zᵣ

		sum(Xᵣ)
	]

  return rtn
end



function Jroessler(V)
	local m = length(V-4)÷6
	local ω,X₀,Xᵣ,Xᵢ,Y₀,Yᵣ,Yᵢ,Z₀,Zᵣ,Zᵢ

	ω = V[end]
	V = reshape(V[1:end-1],2m+1,3)
	X₀,Y₀,Z₀ = V[1,1], V[1,2], V[1,3]
	Xᵣ,Yᵣ,Zᵣ = V[2:2+m-1,1], V[2:2+m-1,2], V[2:2+m-1,3]
	Xᵢ,Yᵢ,Zᵢ = V[2+m:end,1], V[2+m:end,2], V[2+m:end,3]

	local D = ω*(1:m)

	# Jacobian of circular convolution of coefficients of real functions in appropriate format.
	function rCCD(V₀,Vᵣ,Vᵢ)
		local I1,I2,Wᵣ,Wᵢ
		I1 = [ mod(i-j, 2m+1)+1 for i in 0:m, j in 0:m ]
		I2 = [ mod(i+j, 2m+1)+1 for i in 0:m, j in 0:m ]
		Wᵣ = [V₀; Vᵣ; Vᵣ[end:-1:1]]
		Wᵢ = [.0; Vᵢ; -Vᵢ[end:-1:1]]

		return [
			(Wᵣ[I1]+Wᵣ[I2])					(-Wᵢ[I1]+Wᵢ[I2])[:,2:end]
			(Wᵢ[I1]+Wᵢ[I2])[2:end,:]		(Wᵣ[I1]-Wᵣ[I2])[2:end,2:end]
		] / (2m+1)
	end

	local AbyX, AbyY, AbyZ, BbyX, BbyY, BbyZ, CbyX, CbyY, CbyZ, LbyX, LbyY, LbyZ

	AbyX = [
		.0					zeros(1,m)			zeros(1,m)
		zeros(m)			zeros(m,m)			diagm(D)
		zeros(m)			diagm(-D)			zeros(m,m)
	]

	AbyY = [
		-1.0				zeros(1,m)			zeros(1,m)
		zeros(m)			-eye(m)				zeros(m,m)
		zeros(m)			zeros(m,m)			-eye(m)
	]

	AbyZ = [
		-1.0				zeros(1,m)			zeros(1,m)
		zeros(m)			-eye(m,m)			zeros(m,m)
		zeros(m)			zeros(m,m)			-eye(m,m)
	]

	BbyX = [
		1.0					zeros(1,m)			zeros(1,m)
		zeros(m)			eye(m)				zeros(m,m)
		zeros(m)			zeros(m,m)			eye(m)
	]

	BbyY = [
		a					zeros(1,m)			zeros(1,m)
		zeros(m)			a*eye(m)			diagm(D)
		zeros(m)			diagm(-D)			a*eye(m)
	]

	BbyZ = zeros(2m+1,2m+1)

	CbyX = rCCD(Z₀,Zᵣ,Zᵢ)

	CbyY = zeros(2m+1, 2m+1)

	CbyZ = rCCD(X₀,Xᵣ,Xᵢ) + [
		-ℵ					zeros(1,m)			zeros(1,m)
		zeros(m)			-ℵ*eye(m,m)			diagm(D)
		zeros(m)			diagm(-D)			-ℵ*eye(m,m)
	]

	LbyX = [ 1.0			ones(1,m)				zeros(1,m) ]
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

	return [
		AbyX AbyY AbyZ Abyω
		BbyX BbyY BbyZ Bbyω
		CbyX CbyY CbyZ Cbyω
		LbyX LbyY LbyZ .0
	]
end



######

# functions
Dimensions = 3
System = roessler
H = Hroessler
J = Jroessler

# derivative
Epsilon = 1e-2

# cycle finding
TransientIterations = 5000
TransientStepSize = .1
SteadyStateIterations = 25000
SteadyStateStepSize = .01

# initial cycle
InternSamples = 128
ExternSamples = 1024
InterpPrec = 3