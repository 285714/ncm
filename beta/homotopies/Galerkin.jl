#NOTE Galerkin defined this way allows only for numerical derivation.

#TODO dummy plot at Figure()

include("../lib/ncmprojFINDINITIALDATA.jl")
using matsboINTERPOLATE
using matsboNWTN

type Galerkin <: Homotopy
	f::Function
	d::Int

	V::Vector{Float64}

	TIters::Int
	TStepSize::Float64
	SSIters::Int
	SSStepSize::Float64
	m::Int
	Periods::Int

	Galerkin(f::Function, d::Int) = new(f,d)
end

#NOTE untested

function H(G::Galerkin)
	return function H(V::Vector{Float64})
		m = G.m
		ω,ℵ = V[end-2:end]
		C = reduce(hcat, [ complex(V[(i-1)*(2m+1)+(1:m+1)], [0, V[(i-1)*(2m+1)+(m+1)+(1:m)]]) for i in 1:G.d ])
		D = rfft(mapslices(v->G.f([v;ℵ]), irfft(C, 2m, [1]), [2]), [1])
		E = D - im*ω*(0:m) .* C
		return [ reduce(vcat, [ [real(E[:,i]); imag(E[2:end,i])] for i in 1:G.d ]); sum(abs(C[:,1])) ]
	end
end

function J(G::Galerkin)
	error("not defined")
end

global canvasGalerkin #TODO prevent global...

function C(G::Galerkin)
	global canvasGalerkin

	windowGalerkin = @Window("Galerkin Controls", 128, 256, false, true)
	# setproperty!(windowGalerkin, :type_hint, Gtk.GdkWindowTypeHint.TOOLBAR)

	# boxImmerse, toolbarImmerse, canvasImmerse = Immerse.createPlotGuiComponents()
	# delete!(boxImmerse, toolbarImmerse)
	# canvasGalerkin = canvasImmerse

	gridGalerkin = @Grid()
	setproperty!(gridGalerkin, :column_spacing, 5)
	setproperty!(gridGalerkin, :row_spacing, 5)
	push!(windowGalerkin, gridGalerkin)

	canvasGalerkin = @Canvas()
	setproperty!(canvasGalerkin, :expand, true) #necessary for plotting
	gridGalerkin[1,1] = canvasGalerkin

	gridControls = @Grid()
	setproperty!(gridControls, :row_homogeneous, true)
	setproperty!(gridControls, :column_homogeneous, true)
	# setproperty!(gridControls, :vexpand, false)
	# setproperty!(gridControls, :hexpand, true)
	gridGalerkin[1,2] = gridControls

	spinTIters = @SpinButton(0,1e10,1000)
	spinTStepsize = @SpinButton(0,1000,.001)
	spinSSIters = @SpinButton(0,1e10,1000)
	spinSSStepSize = @SpinButton(0,1000,.001)
	spinm = @SpinButton(1,4096,1)
	spinPeriods = @SpinButton(1,128,1)
	buttonFindInitialValue = @Button("Find Initial Value")
	gridControls[1,1] = @Label("Trans. Iterations");	gridControls[2,1] = spinTIters
	gridControls[3,1] = @Label("Trans. StepSize"); 		gridControls[4,1] = spinTStepsize
	gridControls[1,2] = @Label("SS Iterations");		gridControls[2,2] = spinSSIters
	gridControls[3,2] = @Label("SS StepSize"); 			gridControls[4,2] = spinSSStepSize
	gridControls[1,3] = @Label("m");					gridControls[2,3] = spinm
	gridControls[3,3] = @Label("Periods");				gridControls[4,3] = spinPeriods
	gridControls[3:4,4] = buttonFindInitialValue

	# Handlers
	signal_connect(spinTIters, "value_changed") 		do w
		G.TIters = getproperty(w, :value, Int)
	end

	signal_connect(spinTStepsize, "value_changed") 		do w
		G.TStepSize = getproperty(w, :value, Float64)
	end

	signal_connect(spinSSIters, "value_changed") 		do w
		G.SSIters = getproperty(w, :value, Int)
	end

	signal_connect(spinSSStepSize, "value_changed") 	do w
		G.SSStepSize = getproperty(w, :value, Float64)
	end

	signal_connect(spinm, "value_changed")				do w
		G.m = getproperty(w, :value, Float64)
	end

	signal_connect(spinPeriods, "value_changed")		do w
		G.Periods = getproperty(w, :value, Float64)
	end

	signal_connect(buttonFindInitialValue, "clicked") 	do w
		@async handlerFindInitialData(G)
	end

	showall(windowGalerkin)
end



function V(G::Galerkin, V::Vector{Float64})
	local m = ((length(V)-2)÷G.d - 1) ÷ 2
	local C = [ irfft(complex(V[(i-1)*(2m+1)+(1:m+1)], [0, V[(i-1)*(2m+1)+(m+1)+(1:m)]]), 2m) for i in 1:G.d ]

	settings = [
		Geom.line(),
		Theme(background_color=colorant"white"),
		Guide.xticks(ticks=nothing, orientation=:horizontal),
		Guide.yticks(ticks=:auto, orientation=:horizontal),
		Guide.xlabel(""),
		Guide.ylabel(""),
		Coord.cartesian(xmin=.0, xmax=2pi)
	]

	x = linspace(0, 2pi, 2m)
	local P = Plot[ plot(x=x, y=C[i], settings...) for i in 1:G.d ]

	fig = Figure(canvasGalerkin, plot(x=[0], y=[0]))
	fig.cc = vstack(P)

	display(fig)
end



#TODO add missing fields: t0, y0
#TODO add omega
function handlerFindInitialData(G)
	dataT, dataSS, P = findCycle((t,v)->G.f(t,[v;7.0]), .0, rand(G.d), G.TIters, G.TStepSize, G.SSIters, G.SSStepSize)
	cyc,ω = prepareCycle(dataSS, G.SSStepSize, P; fac=G.Periods)

	C = mapslices(cyc, [1]) do v
		tmp = rfft(v)
		m = length(tmp)-1
		f = x -> interpolateTrigonometric(real(tmp[1]), 2*real(tmp[2:end]), -2*imag(tmp[2:end]))(x) / (2m+1)
		tmp = map(f, linspace(.0,2pi,2*G.m+2)[1:end-1])
		tmp = rfft(tmp)
		return [ real(tmp); imag(tmp[2:end]) ]
	end

	C = [reduce(vcat, C); ω; 7.0]
	V(G, C)
	return
	
	# C = newton(H(G), x->forwardDifference(H(G),x), C, predCount(1); callback=(v,h,j)->println(norm(h)))

	#tmp
	# m = G.m
	# C = [ irfft(complex(C[(i-1)*(2m+1)+(1:m+1)], [0, C[(i-1)*(2m+1)+(m+1)+(1:m)]]), 2m) for i in 1:G.d ]

	m = G.m
	ω,ℵ = C[end-2:end]
	C = reduce(hcat, [ complex(C[(i-1)*(2m+1)+(1:m+1)], [0, C[(i-1)*(2m+1)+(m+1)+(1:m)]]) for i in 1:G.d ])
	D = rfft(mapslices(v->G.f([v;ℵ]), irfft(C, 2m, [1]), [2]), [1])
	E = D - im*ω*(0:m) .* C
	F = mapslices(x->irfft(x,2m), E, [1])

	settings = [
		Theme(background_color=colorant"white"),
		Guide.xticks(ticks=nothing, orientation=:horizontal),
		Guide.yticks(ticks=:auto, orientation=:horizontal),
		Guide.xlabel(""),
		Guide.ylabel(""),
		Coord.cartesian(xmin=.0, xmax=2pi)
	]

	x1 = linspace(0, 2pi, size(cyc,1)+1)[1:end-1]
	x2 = linspace(0, 2pi, 2*G.m+1)[1:end-1]
	local Pl = Plot[ plot(
		layer(x=x1, y=cyc[:,i], Geom.line()),
		layer(x=x2, y=F[:,i], Geom.line()), settings...)
		for i in 1:G.d
			]

	fig = Figure(canvasGalerkin, plot(x=[0], y=[0]))
	fig.cc = vstack(Pl)
	display(fig)

	return
	#tmp
end
