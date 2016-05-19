using PyPlot
using Debug
using ODE

include("../src/ncm.jl")
include("../src/diff.jl")

rossler = p -> v -> begin
	a,b,c = p
	x,y,z = v
	return [-y-z ; x + a*y ; b + z*(x-c)]
end


#=
function jacobian(f, ɛ = 0.001)
	return X -> begin
		fX = f(X);
		h  = length(fX)
		w  = length(X)
		J  = zeros(h, w)
		print(" -J:$h×$w- ")
		for i = 1:w
			J[:,i] = ( f(X + ɛ .* ((1:w).==i)) - fX ) ./ ɛ
		end
		J
	end
end
=#



#=
type rossler
	m :: Int			# dimension
	a :: Real
	b :: Real
	c :: Real
	varparam :: Int		# varying parameter (a, b, or c)
	ω :: Real
	cx :: Array{}
	cy :: Array{}
	cz :: Array{}
end

function unwrap(r :: rossler)
	[#=[r.a r.b r.c][r.varparam] ; =#r.ω ; real([r.cx;r.cy;r.cz]) ; imag([r.cx;r.cy;r.cz])]
end

function wrap!(r :: rossler, V)
	#=r.varparam == 1 ? r.a = V[1] : r.varparam == 2 ? r.b = V[1] : r.c = V[1]=#
	V = [0; V] #
	r.ω = V[2]
	k = 2*r.m + 1
	cxr, cyr, czr, cxi, cyi, czi = map(i -> V[i*k+3:(i+1)*k+2], 0:5)
	r.cx = complex(cxr, cxi)
	r.cy = complex(cyr, cyi)
	r.cz = complex(czr, czi)
	r
end

R = rossler(10, 0.1, 0.1, 4, 3, 6.0, rand(21), rand(21), rand(21))


rosslerH = @wrapped R function(r :: rossler)
	k         = -r.m:r.m
	sumMatrix = reverse(r.cx) * transpose(r.cz)
	s         = map(d -> sum(diag(sumMatrix, d)), k)

	eq1 = r.ω .* r.cx .* im .* k  +  r.cy + r.cz							# (1)
	eq2 = r.ω .* r.cy .* im .* k  -  r.cx - r.a .* r.cy						# (2)
	eq3 = r.ω .* r.cz .* im .* k  -  r.b .* (k.==0) + r.c .* r.cz  -  s 	# (3)

	[real([eq1; eq2; eq3]); imag([eq1; eq2; eq3]); sum(real(r.cx))]
	# [real(eq1); imag(eq1); real(eq2); imag(eq2); real(eq3); imag(eq3); sum(real(r.cx))]
end

rosslerDH = jacobian(rosslerH)
=#



include("../mats/lib/matsboUTIL.jl")
include("../mats/lib/matsboNWTN.jl")
include("../mats/lib/matsboRK.jl")
include("../mats/lib/matsboPRED.jl")

include("../mats/lib/ncmprojFINDINITIALDATA.jl")
include("../mats/lib/ncmprojPREPAREDATA.jl")
include("../mats/lib/ncmprojGENERATEFUNCTIONS.jl")

include("../mats/systems.jl")
include("../mats/globals.jl")


c = 8.5
Samples = 256

H, J, rToC, cToR = generateFunctions(System, Dimensions, Epsilon)
Trace, P = findCycle(System, .0, 10*rand(3), TransientIterations, TransientStepSize, SteadyStateIterations, SteadyStateStepSize)
y,ω = prepareCycle(Trace, P, Samples, fac=2)

Y = rfft(y, [1])
V₀ = cToR(Any[matsboUTIL.columns(Y); ω])
V = matsboNWTN.newton(H, matsboNWTN.broyden(H, J), V₀, matsboPRED.predCount(16))


type galerkin
	f :: Function
	param :: Array
	varparam :: Int
	dims :: Int
	N :: Int
	C :: Array
	ω :: Real
end

G = galerkin(rossler, [a,b,c], 3, 3, Samples, [], 0)


function unwrap(g :: galerkin)
	[g.param[g.varparam]; cToR(Any[matsboUTIL.columns(g.C); g.ω])]
end

function wrap!(g :: galerkin, V)
	g.param[g.varparam] = V[1]
	X = rToC(V[2:end]) 
	g.C = first(X)
	g.ω = last(X)[1]
	g
end

galerkinH = @wrapped G function(g :: galerkin)
	local tmp = rfft(mapslices(g.f(g.param), irfft(g.C, g.N, [1]), [2]), [1])
	tmp = broadcast((a,b,c)->a-b*c, tmp, im*g.ω.*(0:g.N÷2), g.C)

	return cToR([ matsboUTIL.columns(tmp); reduce((X,x)->X+real(x), .0, g.C[:,1]) ])
end

galerkinDH = forwardDiff(galerkinH)

wrap!(G, [c; V])


function plotGalerkin(g :: galerkin)
	c = irfft(G.C, G.N, [1])
	c = [c; c[1,:]]
	plot3D(c[:,1], c[:,2], c[:,3])
end

function findLocalMax(g :: galerkin, dim, approx=[])
	if approx == []
		y = irfft(g.C[:,dim], g.N)
		desc = ([y[end]; y[1:end-1]] - y) .< 0
		desc = [desc; desc[1]]
		approx = y[map(i -> desc[i] && !desc[i+1], 1:length(y))]
	end
	approx
end

function timeShift!(g :: galerkin, t)
	g.C .*= exp(-im .* G.ω .* (0:G.N÷2) .* t)
	g
end

function period(g :: galerkin)
	length(findLocalMax(g, 1))
end

function plotIntGalerkin(g :: galerkin, dy = zeros(3), t=4π/g.ω)
	f = (t,x) -> g.f(g.param)(x)
	y0 = irfft(g.C, g.N, [1])[1,:]
	(t, pos) = ode45(f, map(i -> y0[i], 1:length(y0)) + dy, 0:0.001:t)
	x = map(pos -> pos[1], pos)
	y = map(pos -> pos[2], pos)
	z = map(pos -> pos[3], pos)
	subplot(313, projection="3d")
	plot3D(x,y,z)
	p = [pos[end] pos[end÷2]]'
	scatter3D(p[:,1], p[:,2], p[:,3])
	pos
end

function isStable(g :: galerkin, ɛ = 1)
	f = (t,x) -> g.f(g.param)(x)
	y0 = irfft(g.C, g.N, [1])[1,:]
	(t, pos) = ode45(f, map(i -> y0[i], 1:length(y0)), 0:0.001:2π/g.ω)
	norm(pos[end] - pos[1]) < ɛ
end


pc = PredictorCorrector(1, 0.001, 1.8, 2.0, 0.5, 1, (), 0)  #c=6
# pc = PredictorCorrector(1, 0.001, 1.8, 2.0, 0.5, (), 0)  #c=4

#=cm(pc,
   galerkinH,
   galerkinDH,
   m -> (@wrapped G g -> begin
			print("\na=$(g.param[1]), b=$(g.param[2]), c=$(g.param[3]) → ")
			subplot(211, projection="3d")
			plotGalerkin(g)
			x = g.param[3]
			y = g.param[2]
			z = findLocalMax(g, 1)
			subplot(212, projection="3d")
			scatter3D(x .* ones(z), y .* ones(z), z)
			return m.n < 10
		end)(m.v),
   G)=#



t = Task(() -> cm(pc, galerkinH, galerkinDH, produce, G))

maxPlots = Dict()

function runN(n :: Int)
	i = 0
	println("$n steps with ɛ=$(pc.ɛ), κ=$(pc.κ), δ=$(pc.δ), α=$(pc.α):")

	for m in t
		(@wrapped G g -> begin
			print("\n[ a=$(g.param[1]), b=$(g.param[2]), c=$(g.param[3]), ω=$(g.ω) $(isStable(g) ? "✓ " : "")] ")
			subplot(311, projection="3d")
			plotGalerkin(g)

			global maxPlots
			x = g.param[3] # c
			y = g.param[1] # a
			Z = sort(findLocalMax(g, 1))
			p = length(Z)
			merge!(maxPlots, Dict(map(i ->
				("$p.$i", [get(maxPlots,"$p.$i",zeros(0,3)); x y Z[i]]), 1:p)))


			# subplot(212, projection="3d")
			subplot(312)
			# scatter3D(maxPlots[:,1],maxPlots[:,2],maxPlots[:,3])
			for line = values(maxPlots)
				plot(line[:,1], line[:,3])
			end
			xlabel("c")
			ylabel("x-max")
			plotIntGalerkin(g, [0; 1; 0], 4π/g.ω)
		end)(m.v)

		if (i+=1) >= n
			break
		end
	end
end





