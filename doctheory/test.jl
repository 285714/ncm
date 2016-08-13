using PyPlot
ion()


function Jnum(Y, ϵ=1e-4*(1+im))
  Y = reshape(Y,N,n)
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
	tmp = ifft(Y, [1])
	tmp = mapslices(f, tmp, [2])
	tmp = fft(tmp, [1])
	return vec(tmp - im .* ifftshift(-m:m) .* Y)
end

@noinline function J(Y)
  Y = reshape(Y,N,n)
	tmp = ifft(Y, [1])
  M = Array{Complex{Float64}}(N,n,n)
  for i in 1:N
    global asdf = Y[i,:]
    M[i,:,:] = f′(Y[i,:]')
	end
  fft!(M,[1])
  K = im*diagm(ifftshift(-m:m))
  tmp = [ circulant(M[:,k,l]./N) - (k==l ? K : 0) for k in 1:n, l in 1:n ]
  tmp = reducedim(hcat, tmp, [2], Array{Complex}(N,0))
  tmp = reducedim(vcat, tmp, [1], Array{Complex}(0,n*N))[1]
end

maxabs(Jnum(Y) - J(Y))
matshow(real(Jnum(Y) - J(Y))); colorbar()
matshow(imag(Jnum(Y) - J(Y))); colorbar()
matshow(real(J(Y))); colorbar()
matshow(real(Jnum(Y))); colorbar()
close()

function circulant(v)
  N = length(v)
  return [ v[mod(i-j,N)+1] for i in 1:N, j in 1:N ]
end

circulant(rand(3))
circulant([1,2,3])

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

f(x) = lorenz(0,[x;28])
f′(x) = lorenz′(0,[x;28])[:,1:end-1]


m = 31
N = 2m+1
n = 3
Y = 100*rand(Complex{Float64}, N*n)

H(Y)
J(Y)
Jnum(Y)


newton(f,f′,x) = x - inv(f′(x))*f(x)

Y = 5*rand(Complex{Float64}, N*n)
Ytmp = Y
Ytmp = newton(H,J,Ytmp); println(norm(H(Ytmp)))
Ytmp = newton(H,Jnum,Ytmp); println(norm(H(Ytmp)))
H(Ytmp)
J(Ytmp)
