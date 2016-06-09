#NOTE Galerkin defined this way allows only for numerical derivation.

#TODO dummy plot at Figure()

using Immerse
using Gtk.ShortNames
include("../types.jl")

type Galerkin <: Homotopy
	f::Function
	d::Int
end

#NOTE untested

function H(G::Galerkin)
	local m = ((length(V)-2)÷G.d - 1) ÷ 2
	return function H(V::Vector{Float64})
		local C = reduce(hcat, [ complex(V[(i-1)*(2m+1)+(1:m+1)], [0, V[(i-1)*(2m+1)+(m+1)+(1:m)]]) for i in 1:G.d ])
		C = rfft(mapslices(G.f, irfft(C, 2m, [1]), [2]), [1])
		return reduce(vcat, [ [real(C[i]); imag(C[i][2:end])] for i in 1:G.d ])
	end
end

function J(G::Galerkin)
	error("not defined")
end

#TODO labels, sizes, plot

global canvasGalerkin

function C(G::Galerkin)
	global canvasGalerkin

	windowGalerkin = @Window("Galerkin Controls", 128, 256, false, true)
	# setproperty!(windowGalerkin, :type_hint, Gtk.GdkWindowTypeHint.TOOLBAR)

	# boxImmerse, toolbarImmerse, canvasImmerse = Immerse.createPlotGuiComponents()
	# delete!(boxImmerse, toolbarImmerse)
	# canvasGalerkin = canvasImmerse

	gridGalerkin = @Grid()
	push!(windowGalerkin, gridGalerkin)

	canvasGalerkin = @Canvas()
	setproperty!(canvasGalerkin, :expand, true) #necessary for plotting
	gridGalerkin[1,1] = canvasGalerkin

	gridControls = @Grid()
	setproperty!(gridControls, :vexpand, false)
	setproperty!(gridControls, :hexpand, true)
	setproperty!(gridControls, :row_homogeneous, true)
	setproperty!(gridControls, :column_homogeneous, true)
	gridGalerkin[1,2] = gridControls

	spinTIters = @SpinButton(0,1e10,1000)
	spinTStepsize = @SpinButton(0,1000,.001)
	spinSSIters = @SpinButton(0,1e10,1000)
	spinSSStepSize = @SpinButton(0,1000,.001)
	buttonFindInitialValue = @Button("Find Initial Value")
	gridControls[1,1] = spinTIters
	gridControls[2,1] = spinTStepsize
	gridControls[1,2] = spinSSIters
	gridControls[2,2] = spinSSStepSize
	gridControls[1:2,3] = buttonFindInitialValue

	buttonDoublePeriod = @Button("Double Period")
	buttonHalfPeriod = @Button("Half Period")
	gridControls[4,1] = buttonDoublePeriod
	gridControls[4,2] = buttonHalfPeriod

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
