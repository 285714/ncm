#simplified concept
#	two components: continuation & system
#	no interface requirements: components call stuff and it either exists or it doesn't
#		exceptions: continuationExec, systemExec must exist
#		types only for datatypes (if at all), generalization is user responsibility
#		global Proj variable is accessed by both components
#			never simultaneously, no sync, no communication (in general, otherwise: user responsibility)
#		no contraints on what components are allowed to do: user responsibility

#	this file specifies the whole interface...

#TODO ifft for trigoninterp

addprocs(1)

push!(LOAD_PATH, "$(pwd())/lib")
using Gtk.ShortNames
@everywhere using mbNewton, mbPred, mbUtil, mbInterpolate, mbObserve
@everywhere include("Model.jl")
include("lib/ncmprojMKCONTROLGRID.jl")

#(do menu stuff, other global GUI stuff, saving, ...)

#create empty project: vector of branches of solutions
global project = Project()
global lockProject = ReentrantLock()

save(filename) = open("save/$(filename)", "w") do f serialize(f, project) end
load(filename) = (global project = open(deserialize, "save/$(filename)"))

#select continuation method, select system
# include(open_dialog("Select continuation method.", Gtk.GtkNullContainer(), ASCIIString[]))
# include(open_dialog("Select system.", Gtk.GtkNullContainer(), ASCIIString[]))

# @everywhere include("system/lorenz.jl")
@everywhere include("system/roessler.jl")

include("lib/ncmprojBIFPLOT.jl") #needs thread 1
@schedule B = BifPlot(project)

@everywhere include("continuation/PC.jl")
@schedule pcGUI()

@everywhere include("lib/ncmprojGALERKIN.jl")
@schedule galerkinGUI()


#TODO-2DAY:
#	abstraction: plotting ✓, galerkin, continuation, project (! - for access, sync, sanity, ...)
#		lock
#		project: all operations
#	howto associate additional info with solutions/branches? transparent, general...
#		projects must be method agnostic...
#		plots, projection, stepsize, ...
#		encapsulate to simulate inheritance? serializability?

# (perturbation)
# bug: stepsize is associated with branch, not correct on deletion... reset, make changeable, associate with solution
# bug: apparently activeSoltion read from wrong (global) project...

# poincare ✓
# improv RK	?
# perturbation ✓
# lorenz

#TODO keys: delete from here, delete till here
#TODO h box, fix
#TODO broyden
#TODO mark current branch


# plot:
#	mark curr branch ✓
#	color ✓
#	lock branch (make unselectable, how to unlock?, how switch branches efficiently?)
#	rebuild completely as observable


# winkel pred ✓
# speichern -> encapsulate galerkin, make serializable
# lorenz
# roessler bif

#mark perturb branches
#auto new branch
#singleton perturb branch?
#new panel?
#special perturb mode?
#change perturb style... make factor for rand

#broyden pc, keep J per branch


##=
include("lib/ncmprojGALERKINHELPER.jl")

c = 6.0
m = 16
cyc, P, ω = findCyclePoincare((t,v)->f(t,[v;c]), rand(3),
	nIntersections=120, maxCycles=30, sampleSize=m,
	transientIterations=4000, transientStepSize=.1,
	steadyStateStepSize=.1)

C = resample(rfft(cyc, [1]), m)
C = [vec(vcat(real(C), imag(C[2:end,:]))); ω]

D = [C;c]
setActiveSolution(project, D)
@time D = newton(H, J2, D, predCount(100), callback=(v,H,J)->println(norm(v)))

J2 = f′ToJ(f′)
J4(D) = forwardDifference(H,D)

ion()
D = 1000*rand(107)
matshow(J(D)-J2(D)); colorbar()
matshow(J4(D)-J2(D)); colorbar()
matshow(J2(D)); colorbar()
matshow(J2(D)/J(D)); colorbar()
matshow(map(x->clamp(x,-10,10),J(D))); colorbar()
matshow(map(x->clamp(x,-10,10),J2(D))); colorbar()
matshow(map(x->clamp(x,-100,100),J2(D)-J(D))[2*(2m+1)+(1:2m+1),2*(2m+1)+(1:2m+1)]); colorbar()
matshow(map(x->clamp(x,-10,10),J(D))[2*(2m+1)+(1:2m+1),2*(2m+1)+(1:2m+1)]); colorbar()
matshow(map(x->clamp(x,-10,10),J2(D))[2*(2m+1)+(1:2m+1),2*(2m+1)+(1:2m+1)]); colorbar()
matshow(J2(D)); colorbar()


tmp = J(D)
tmp = tmp[2*(2m+1)+(1:2m+1),2*(2m+1)+(1:2m+1)]
a = tmp[m+2:end,1:m+1]
b = tmp[1:m+1,m+2:end]

a,b = .1,.1
=#



#######

V = D

m = length(V-5)÷6
ω,ρ = V[end-1:end]
V = reshape(V[1:end-2],2m+1,3)
V = complex(V[1:m+1,:], [zeros(3)'; V[m+2:end,:]])

T = irfft(V, 2m, [1])
M = Array{Float64}(size(T,1), size(T,2), size(T,2)+1) #n, comp, deriv
for i in 1:size(T,1) M[i,:,:] = f′([T[i,:]'; ρ]) end

shiftMatrix(V) = [ V[mod(i-j,length(V))+1] for i in 1:length(V), j in 1:length(V) ]
M′ = [ shiftMatrix(rfft(M[:,i,j]/(2m))) for i in 1:size(T,2), j in 1:size(T,2) ]
M′ = [ rCCD(rfft(M[:,i,j])) for i in 1:size(T,2), j in 1:size(T,2) ]

for i in 1:size(T,2)
	D = diagm(ω*(0:m))
	M′[i,i][m+2:end,1:m+1] -= D[2:end,:]
	M′[i,i][1:m+1, m+2:end] += D[:,2:end]
end

for i in 1:size(T,2), j in 1:size(T,2)
	A,B = real(M′[i,j]), imag(M′[i,j])
	# i==j && (B[diagind(B)] -= ω*(0:m))
	M′[i,j] = [
		A			B'[:,2:end]
		B[2:end,:]	A[2:end,2:end]
	]
end

M′ = reducedim(vcat, M′, [1], Array{Float64}(0,size(M′[1,1],2)))
M′ = reducedim(hcat, M′, [2], Array{Float64}(size(M′[1,1],1),0))[1]

ddω = -im*(0:m).*V
ddω = vec([real(ddω); imag(ddω)[2:end,:]])

ddρ = reduce(hcat, [ rfft(M[:,i,size(T,2)+1]) for i in 1:size(T,2) ])
ddρ = vec([real(ddρ); imag(ddρ)[2:end,:]])

return [
	[M′ ddω ddρ]
	[ones(1,m+1) zeros(1,size(M′,2)-m+1)]
	]




#####


V = D

m = length(V-5)÷6

ω,ℵ = V[end-1:end]
V = reshape(V[1:end-2],2m+1,3)
X₀,Y₀,Z₀ = V[1,1], V[1,2], V[1,3]
Xᵣ,Yᵣ,Zᵣ = V[2:2+m-1,1], V[2:2+m-1,2], V[2:2+m-1,3]
Xᵢ,Yᵢ,Zᵢ = V[2+m:end,1], V[2+m:end,2], V[2+m:end,3]

D′ = ω*(1:m)

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

AbyX = [
	.0					zeros(1,m)			zeros(1,m)
	zeros(m)			zeros(m,m)			diagm(D′)
	zeros(m)			diagm(-D′)			zeros(m,m)
	]

AbyY = AbyZ = -eye(2m+1)

BbyX = eye(2m+1)

BbyY = [
	a					zeros(1,m)			zeros(1,m)
	zeros(m)			a*eye(m)			diagm(D′)
	zeros(m)			diagm(-D′)			a*eye(m)
	]

BbyZ = zeros(2m+1,2m+1)

CbyX = rCCD(Z₀,Zᵣ,Zᵢ)

CbyY = zeros(2m+1, 2m+1)

CbyZ = rCCD(X₀,Xᵣ,Xᵢ) + [
	-ℵ					zeros(1,m)			zeros(1,m)
	zeros(m)			-ℵ*eye(m,m)			diagm(D′)
	zeros(m)			diagm(-D′)			-ℵ*eye(m,m)
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

A = [
	AbyX AbyY AbyZ Abyω Abyℵ
	BbyX BbyY BbyZ Bbyω Bbyℵ
	CbyX CbyY CbyZ Cbyω Cbyℵ
	LbyX LbyY LbyZ Lbyω Lbyℵ
	]




test = rCCD(X₀,Xᵣ,Xᵢ) + [
	-ℵ					zeros(1,m)			zeros(1,m)
	zeros(m)			-ℵ*eye(m,m)			zeros(m,m)
	zeros(m)			zeros(m,m)			-ℵ*eye(m,m)
	]
test = CbyZ

matshow(test); colorbar()
matshow(M′[3,3]); colorbar()
matshow(M′[3,3]-test); colorbar()
