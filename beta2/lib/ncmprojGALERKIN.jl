using Gtk.ShortNames
using matsboINTERPOLATE, matsboNWTN, matsboPRED, matsboUTIL
include("$(pwd())/lib/ncmprojFINDINITIALDATA.jl")
include("$(pwd())/lib/ncmprojMKCONTROLGRID.jl")

function galerkinGUI()
	windowGalerkin = @Window("Galerkin Controls", 256, 256, false, true)
	setproperty!(windowGalerkin, :resizable, false)
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
	signal_connect(w -> @async(handlerFindInitialData(w, dataItSS)), buttonFindInitialValue, "clicked")
	gridGalerkin[1,2] = buttonFindInitialValue

	showall(windowGalerkin)
end



function handlerFindInitialData(w, dataGUI)
	setproperty!(w, :sensitive, false)
	#TODO function/macro bringIntoScope(D::Dict)
	TIters, SSIters, TStepSize, SSStepSize, Periods, m, c₀ = map(x->dataGUI[x], ["Trans. Iterations", "SS Iterations", "Trans. StepSize", "SS StepSize", "Periods", "m", "c₀"])

	C = @fetch begin
		dataT, dataSS, P = findCycle((t,v)->f(t,[v;c₀]), .0, rand(3), TIters, TStepSize, SSIters, SSStepSize)
		cyc,ω = prepareCycle(dataSS, SSStepSize, P; fac=Periods)

		local C = mapslices(cyc, [1]) do v
			tmp = rfft(v)
			tmp = map(linspace(.0,2pi,2*m+2)[1:end-1]) do x
				interpolateTrigonometric(real(tmp[1]), 2*real(tmp[2:end]), -2*imag(tmp[2:end]))(x) / (2m+1)
			end
			tmp = rfft(tmp)
			return [ real(tmp); imag(tmp[2:end]) ]
		end

		C = [reduce(vcat, C); ω]
		Htmp(V) = H([V; c₀]) # R^N -> R^N system
		Jtmp(V) = J([V; c₀])[:, 1:end-1] # R^N -> R^(NxN) system
		newton(Htmp, Jtmp, C, predCount(10) ∧ predEps(1e-10))
	end

	global project
	project.activeSolution = [C; c₀]

	setproperty!(w, :sensitive, true)
	return Void
end



function toTrajInterp(V, d)
	m = (length(V)-d-2)÷(2*d)
	tmp = reshape(V[1:end-2], 2m+1, d)
	return matsboUTIL.vectorize(x -> [ matsboINTERPOLATE.interpolateTrigonometric(tmp[1,i], 2tmp[2:2+m-1,i], -2tmp[2+m:end,i])(x) / (2m+1) for i in 1:d ])
end


function projection(V)
	f = toTrajInterp(V,3)
	rtn = Float64[]

	dt = 2pi/1025 #TODO ...
	for t in .0:dt:2pi
		if f(t)[1] ≤ .0 ≤ f(t+dt)[1]
			x = matsboUTIL.bisection(x->f(x)[1], t, t+dt, ϵ=1e-4) #TODO precision
			push!(rtn, norm(f(x)))
		end
	end

	return rtn
end
