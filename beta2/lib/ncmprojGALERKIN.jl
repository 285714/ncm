type Galerkin
	f::Function
	H::Function
	J::Function
end

using Gtk.ShortNames
using matsboINTERPOLATE, matsboNWTN, matsboPRED, matsboUTIL

include("../lib/ncmprojFINDINITIALDATA.jl")
include("../lib/ncmprojMKCONTROLGRID.jl")

function systemExec()
	roesslerGUI()
end

function galerkinGUI()
	windowGalerkin = @Window("Galerkin Controls", 128, 256, false, true)
	# setproperty!(windowGalerkin, :type_hint, Gtk.GdkWindowTypeHint.TOOLBAR)

	gridGalerkin = @Grid()
	setproperty!(gridGalerkin, :column_spacing, 5)
	setproperty!(gridGalerkin, :row_spacing, 5)
	push!(windowGalerkin, gridGalerkin)

	local dataItSS, gridItSS
	local initItSS = [
		("Trans. Iterations", Int, 0, 1e8, 1000),
		("SS Iterations", Int, 0, 1e8, 1000),
		("Trans. StepSize", Float64, .0, 1.0, .001),
		("SS StepSize", Float64, 0, 1.0, .001),
		("c₀", Float64, .0, 100.0, .1),
		("Periods", Int, 1, 128,1),
		("m", Int, 0, 4096, 1)
	]
	dataItSS, gridItSS = mkControlGrid(initItSS, 1)
	gridGalerkin[1,1] = gridItSS

	buttonFindInitialValue = @Button("Find Initial Solution")
	signal_connect(w -> @async(handlerFindInitialData(dataItSS)), buttonFindInitialValue, "clicked")
	gridGalerkin[1,2] = buttonFindInitialValue

	showall(windowGalerkin)
end



function handlerFindInitialData(dataGUI)
	#TODO function/macro bringIntoScope(D::Dict)
	TIters, SSIters, TStepSize, SSStepSize, Periods, m, c₀ = map(x->dataGUI[x], ["Trans. Iterations", "SS Iterations", "Trans. StepSize", "SS StepSize", "Periods", "m", "c₀"])

	dataT, dataSS, P = findCycle((t,v)->roessler(t,[v;c₀]), .0, rand(3), TIters, TStepSize, SSIters, SSStepSize)
	cyc,ω = prepareCycle(dataSS, SSStepSize, P; fac=Periods)

	local C = mapslices(cyc, [1]) do v
		tmp = rfft(v)
		tmp = map(linspace(.0,2pi,2*m+2)[1:end-1]) do x
			interpolateTrigonometric(real(tmp[1]), 2*real(tmp[2:end]), -2*imag(tmp[2:end]))(x) / (2m+1)
		end
		tmp = rfft(tmp)
		return [ real(tmp); imag(tmp[2:end]) ]
	end

	C = [reduce(vcat, C); ω; c₀]
	C = newton(H, J, C, predCount(10) ∧ predEps(1e-10))

	global project
	B = Branch(C)
	B.D["hUp"], B.D["hDown"] = 1.0, -1.0
	push!(project.branches, B)
end



function toTrajInterp(V, d)
	m = (length(V)-d-2)÷(2*d)
	tmp = reshape(V[1:end-2], 2m+1, d)
	return matsboUTIL.vectorize(x -> [ matsboINTERPOLATE.interpolateTrigonometric(tmp[1,i], 2tmp[2:2+m-1,i], -2tmp[2+m:end,i])(x) / (2m+1) for i in 1:d ])
end


function galerkinProjection(V)
	f = toTrajInterp(V,3)
	rtn = Float64[]

	dt = 2pi/1025
	for t in .0:dt:2pi
		if f(t)[1] ≤ .0 ≤ f(t+dt)[1]
			x = matsboUTIL.bisection(x->f(x)[1], t, t+dt) #TODO precision
			push!(rtn, norm(f(x)))
		end
	end

	return rtn
end
