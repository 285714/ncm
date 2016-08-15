f(x) = [
  -x[2]-x[3]
  x[1] + .1*x[2]
  .1 + x[3]*(x[1]-14)
  ]

f′(x) = [
  0     -1      -1
  1     .1      0
  x[3]   0      (x[1]-14)
  ]

m = 31
N = 2m+1
n = 3
Y = 1000*rand(Complex{Float64}, N*n)

using PyPlot
ion()

function circulant(v)
  N = length(v)
  return [ v[mod(i-j,N)+1] for i in 1:N, j in 1:N ]
end

function Jnum(Y, ϵ=1e-4)
	h₀ = H(Y)
	tmp = cell(length(Y))
	for i in 1:length(Y)
		Y[i] += ϵ
		tmp[i] = vec((H(Y) - h₀) ./ ϵ)
		Y[i] -= ϵ
	end
	return reduce(hcat, tmp)
end

function H(Y)
	Y = reshape(Y,N,n)
	# Y = [Y; conj(Y[end:-1:2,:])]
	tmp = ifft(Y, [1])
	tmp = mapslices(f, tmp, [2])
	tmp = fft(tmp, [1])
	# return vec(tmp)
	return vec(tmp - im .* ifftshift(-m:m) .* Y)
end

@noinline function J(Y)
	Y = reshape(Y,N,n)
	tmp = ifft(Y, [1])
	M = Array{Complex{Float64}}(N,n,n)
	for i in 1:N
		M[i,:,:] = f′(tmp[i,:]')
	end
	ifft!(M,[1])
	K = spdiagm(im*ifftshift(-m:m))
	tmp = [ circulant(M[:,k,l]) - (k==l ? K : 0) for k in 1:n, l in 1:n ]
	tmp = reducedim(hcat, tmp, [2], Array{Complex}(N,0))
	tmp = reducedim(vcat, tmp, [1], Array{Complex}(0,n*N))[1]
end

matshow(real(Jnum(Y) - J(Y))); colorbar()
matshow(imag(Jnum(Y) - J(Y))); colorbar()
matshow(abs(Jnum(Y)) ./ (abs(J(Y)))+1); colorbar()
matshow(real(J(Y))); colorbar()
matshow(real(Jnum(Y))); colorbar()
matshow(imag(J(Y))); colorbar()
matshow(imag(Jnum(Y))); colorbar()

close()

newton(f,f′,x) = x - inv(f′(x))*f(x)

Ytmpa = 5*rand(Complex{Float64}, N*n)
Ytmp = Ytmpa
Ytmp = newton(H,J,Ytmp); println(norm(H(Ytmp)))
Ytmp = newton(H,Jnum,Ytmp); println(norm(H(Ytmp)))
H(Ytmp)
J(Ytmp)

@noinline function J2(Y)
	Y = reshape(Y,N,n)
	tmp = ifft(Y, [1])
	M = Array{Complex{Float64}}(N,n,n)
	for i in 1:N
		M[i,:,:] = f′(tmp[i,:]')
	end
	ifft!(M,[1])
	K = spdiagm(im*ifftshift(-m:m))
	tmp = [ circulant(M[:,k,l]) - (k==l ? K : 0) for k in 1:n, l in 1:n ]
end

function realJ(Y)
	Y = reshape(Y,m+1+m,n)
	Y = Y[1:m+1,:] + im*[zeros(1,n); Y[m+2:end,:]]
	Y = [Y; conj(Y[end:-1:2,:])]
	tmp = J2(Y)
	tmp = map(x->x[1:m+1,1:m+1], tmp)
	tmp = map(x->[real(x) -imag(x)[:,2:end]; imag(x)[2:end,:] real(x)[2:end,2:end]], tmp)
	tmp = reducedim(hcat, tmp, [2], Array{Float64}(N,0))
	tmp = reducedim(vcat, tmp, [1], Array{Float64}(0,N*n))[1]
end

Y = 1000000000*rand(N*n)
tmp1 = realJ(Y)
tmp2 = Jgen(Y)
tmp3 = Jroessler(Y)

matshow(tmp2-tmp3); colorbar()
matshow(abs(tmp1)./(abs(tmp3)+.01)); colorbar()
matshow(tmp3); colorbar()




function Jroessler(V)
	local α=.1, β=.1, ℵ=14
	local X₀,Xᵣ,Xᵢ,Y₀,Yᵣ,Yᵢ,Z₀,Zᵣ,Zᵢ
	V = reshape(V,N,n)

	X₀,Y₀,Z₀ = V[1,1], V[1,2], V[1,3]
	Xᵣ,Yᵣ,Zᵣ = V[2:2+m-1,1], V[2:2+m-1,2], V[2:2+m-1,3]
	Xᵢ,Yᵢ,Zᵢ = V[2+m:end,1], V[2+m:end,2], V[2+m:end,3]

	local D = (1:m)

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
		AbyX AbyY AbyZ
		BbyX BbyY BbyZ
		CbyX CbyY CbyZ
	]
end





function fToH(f)
	return function H(V::Vector{Float64})
		V = reshape(V,N,n)
		anchor = sum(V[1:m+1]) #TODO hm...

		ω,ρ = V[end-1:end]
		V = reshape(V[1:end-2],2m+1,3)
		V = complex(V[1:m+1,:], [zeros(3)'; V[m+2:end,:]])

		defect = rfft(mapslices(x->f(0,[x;ρ]), irfft(V, 2m, [1]), [2]), [1]) - im*ω*(0:m) .* V
		defect = vec([ real(defect); imag(defect)[2:end,:] ])

		return [defect; anchor]
	end
end

function f′ToJ(f′)
	return function J(V::Vector{Float64})
		V = reshape(V,N,n)
		V = complex(V[1:m+1,:], [zeros(1,n);V[m+2:end,:]])
		T = irfft(V, N, [1])

		M = Array{Float64}(size(T,1), size(T,2), size(T,2)) #n, comp, deriv
		for i in 1:N
			M[i,:,:] = f′(T[i,:]')
		end

		M′ = [ rCCD(rfft(M[:,i,j])) for i in 1:size(T,2), j in 1:size(T,2) ]

		for i in 1:size(T,2)
			D = diagm(0:m)
			M′[i,i][m+2:end,1:m+1] -= D[2:end,:]
			M′[i,i][1:m+1, m+2:end] += D[:,2:end]
		end

		M′ = reducedim(vcat, M′, [1], Array{Float64}(0,size(M′[1,1],2)))
		M′ = reducedim(hcat, M′, [2], Array{Float64}(size(M′[1,1],1),0))[1]

		return M′
	end
end

Jgen = f′ToJ(f′)

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
