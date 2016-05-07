using PyPlot
include("lib/matsboUTIL.jl")
include("lib/matsboNWTN.jl")
include("lib/matsboRK.jl")
include("lib/matsboPRED.jl")

include("lib/ncmprojFINDINITIALDATA.jl")
include("lib/ncmprojHGENERIC.jl")

include("systems.jl")

#	convert input/output vector
function rToC(V)
  local N = (length(V)-1)÷3
  local c(x) = [ x[1]+0im; complex(x[2:end÷2+1], x[end÷2+2:end]) ]
  return matsboUTIL.inputConvert(V, [1;N+1;2N+1;3N+1;3N+2], [c,c,c,identity])
end

function cToR(V...)
  local d(x) = [real(x[1]); real(x[2:end]); imag(x[2:end])]
  return matsboUTIL.outputConvert([d,d,d,identity], V...)
end

# optimization target for roessler
#		operating on real Fourier coefficients plus ω
function Hroessler(V::Vector{Float64})
  local X,Y,Z,ω
  X,Y,Z,ω = rToC(V)
  local N = (length(X)-1)*2

  local A,B,C
	local D = im*ω.*(0:N÷2) #deriv
  A = -Y - Z - D.*X
  B = X + a*Y - D.*Y
  C = [N*b; zeros(N÷2)] - c*Z + rfft(irfft(Z,N).*irfft(X,N)) - D.*Z

  return cToR(A,B,C,reduce((X,x)->X+real(x), .0, X))
end

# derivative
Jroessler(v) = matsboNWTN.forwardDifference(Hroessler, v, ϵ=1e-4)


# generic objective function idea...
# function Htest(V)
# 	local X,Y,Z,ω
# 	X,Y,Z,ω = rToC(V)
# 	local N = (length(X)-1)*2
#
# 	local tmp = mapslices(roessler, irfft([X Y Z], N, [1]), [2])
# 	local tmp2 = tmp[1,1] # not working?! why???
# 	tmp = rfft(tmp, [1])
#
# 	local A,B,C
# 	local D = im*ω.*(0:N÷2)
# 	A = tmp[:,1]-D.*X
# 	B = tmp[:,2]-D.*Y
# 	C = tmp[:,3]-D.*Z
#
# 	return cToR(A,B,C,tmp2)
# end
#
# Jtest(v) = matsboNWTN.forwardDifference(Htest, v, ϵ=1e-4)


# callback for interactive plotting
global ax1, ax2, distance
function initPlot()
		global ax1, ax2, distance
		distance = Float64[]
    ion(); figure()
		global ax1 = subplot2grid((3,1), (0, 0), rowspan=2, projection="3d"); plot([.0],[.0],[.0])
		global ax2 = subplot2grid((3,1), (2, 0))
end

function callbackPlot(V, H, J)
	global ax1, ax2, distance

  local X,Y,Z,ω,N
  X,Y,Z,ω = rToC(V); N = 2*(length(X)-1)
  sca(ax1); ax1["lines"][end]["set_color"]("gray"); hold(true); plot(irfft(X,N), irfft(Y,N), irfft(Z,N), color="r")

	push!(distance, sumabs(H))
	sca(ax2); hold(false); plot(distance)

	draw(); sleep(.05)
end






### Start ###

# choose c...
#		4.0, 6.0, 8.5, 12.0 working
#		8.7, 12.6 (period doubling) working only if period is doubled...
c = 12.6

# find representative cycle
TransientIterations, TransientStepSize = 5000, .1
SteadyStateIterations, SteadyStateStepSize = 25000, .01
data, P = findRepresentativeCycle(roessler, .0, rand(3), TransientIterations, TransientStepSize, SteadyStateIterations, SteadyStateStepSize)

# cut out singe period
ω = 2pi / (P * SteadyStateStepSize)
cyc = data[end-P:end,:]

# good double period
# 	this appears to work =D
#		TODO think about how to automate
#ω = 2pi / (2P * SteadyStateStepSize)
#cyc = data[end-2P:end,:]

# simple doubling
# ω /= 2
# cyc = [cyc; cyc]

# resample
y = mapslices(V->matsboUTIL.interpolate(V, 3, linspace(1.0, size(cyc,1), 512)), cyc, [1])


# rotate s.t. x(0)≈0
# this (↓) is ultimately important! optimization is super sensitive to it...
# TODO maybe scale last equation in target somehow?
y = circshift(y,[-findmin(abs(y[:,1]))[2]+1])

# plots of approx steady-state trajectory, extracted cycle, resampled/shifted version
#ioff(); figure(); subplot(111,projection="3d"); hold(true)
#plot(data[:,1], data[:,2], data[:,3], color="b")
#plot(cyc[:,1], cyc[:,2], cyc[:,3], color="r")
#plot(y[:,1], y[:,2], y[:,3], color="g")
#show()

# generate starting vector
Y = rfft(y,[1])
V₀ = cToR(Y[:,1], Y[:,2], Y[:,3], ω)

# test
include("lib/matsboUTIL.jl")
include("lib/ncmprojHGENERIC.jl")
Htest = Hgeneric(roessler, 3)
Jtest(v) = matsboNWTN.forwardDifference(Htest, v, ϵ=1e-4)

# optimize
#@time V₁ = matsboNWTN.newton(Hroessler, Jroessler, V₀, matsboPRED.predCount(16); init=initPlot, callback=callbackPlot)
@time V₂ = matsboNWTN.newton(Hroessler, matsboNWTN.broyden(Hroessler, Jroessler), V₀, matsboPRED.predCount(32); init=initPlot, callback=callbackPlot)
#@time V₃ = matsboNWTN.newton(Htest, Jtest, V₀, matsboPRED.predCount(16); init=initPlot, callback=callbackPlot)
@time V₄ = matsboNWTN.newton(Htest, matsboNWTN.broyden(Htest, Jtest), V₀, matsboPRED.predCount(32); init=initPlot, callback=callbackPlot)

sumabs(Hroessler(V₀))
sumabs(Hroessler(V₁))
sumabs(Hroessler(V₂))
sumabs(Hroessler(V₃))
sumabs(Hroessler(V₄))
