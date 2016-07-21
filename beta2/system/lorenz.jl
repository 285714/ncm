# stable periodic trajectories for ρ ∈ {350, 160, 100.5, 9.65}
include("../lib/ncmprojGALERKINHELPER.jl")

global σ=10, β=8/3
lorenz(t,v) = [
	σ*(v[2]-v[1])
	v[1]*(v[4]-v[3])-v[2]
	v[1]*v[2]-β*v[3]
	]
lorenz(v) = lorenz(0, v)

lorenz′(t,v) = [
	-σ σ 0 0
	(v[4]-v[3]) -1 -v[1] v[1]
	v[2] v[1] -β 0
	]
lorenz′(v) = lorenz′(0, v)

f, f′, H, J = lorenz, lorenz′, fToH(lorenz), f′ToJ(lorenz′)



#=
# DFT(N) = [ exp(-im*2pi/N * i*j) for i in 0:N-1, j in 0:N-1 ] ./ sqrt(N)

# f,f′,H,J = roessler, roessler′, Hroessler, Jroessler
# f


using PyPlot, mbNewton

m = 10
W = V = 1*rand(m*3*2+3+2)

m = length(V-5)÷6
anchor = sum(V[1:m+1])
ω,ρ = V[end-1:end]
V = reshape(V[1:end-2],2m+1,3)
V = complex(V[1:m+1,:], [zeros(3)'; V[m+2:end,:]])
T = irfft(V, 2m, [1])
M = Array{Float64}(size(T,1), size(T,2), size(T,2)+1) #n, comp, deriv
for i in 1:size(T,1) M[i,:,:] = f′([T[i,:]'; ρ]) end

shiftMatrix(V) = [ V[mod(i-j,length(V))+1] for i in 1:length(V), j in 1:length(V) ]

# M′ = [ (F*broadcast(*,M[:,i,j],F′))[1:m+1, 1:m+1] for i in 1:size(T,2), j in 1:size(T,2) ] #normal matrices!
M′ = [ shiftMatrix(ifft(M[:,i,j]))[1:m+1, 1:m+1] for i in 1:size(T,2), j in 1:size(T,2) ]
# M′2 = [ fft(transpose(ifft(diagm(M[:,i,j]),[1])),[1])[1:m+1, 1:m+1] for i in 1:size(T,2), j in 1:size(T,2) ] #TODO...

for i in 1:size(T,2), j in 1:size(T,2)
	# M′[i,j][diagind(M′[i,j])[2:end-1]] += M′[i,j][1,1] #TODO why?!
	if i == j
		M′[i,i][diagind(M′[i,i])] -= im*ω*(0:m)
	end
end

map!(M′) do tmp
	A,B = real(tmp), imag(tmp)
	rtn = [
		A			-B[:,2:end]
		B[2:end,:] 	A[2:end,2:end]
	]
	# rtn[end,end] = 0
	return rtn
end

ddρ = [ rfft(M[:,i,size(T,2)+1]) for i in 1:size(T,2) ]
ddρ = reduce(vcat, map(x->[real(x);imag(x)[2:end]], ddρ))


M′ = reducedim(vcat, M′, [1], Array{Float64}(0,size(M′[1,1],2)))
M′ = reducedim(hcat, M′, [2], Array{Float64}(size(M′[1,1],1),0))[1]

G = forwardDifference(H,W)[1:end-1,1:end-2]
# G = J(W)[1:end-1,1:end-2]
matshow(G); colorbar(); matshow(M′); colorbar()
matshow(G - M′); colorbar()
matshow(G / M′); colorbar()




J′ = f′ToJ(f′)
A,B,C = J(V), J′(V), forwardDifference(H,W)

std(A-B)
std(A-C)
std(B-C)

matshow(A)
matshow(B)
matshow(A-B); colorbar()
matshow(A-C); colorbar()
=#
