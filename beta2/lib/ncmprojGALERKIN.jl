using Gtk.ShortNames
using matsboINTERPOLATE, matsboNWTN, matsboPRED, matsboUTIL
include("$(pwd())/lib/ncmprojFINDINITIALDATA.jl")
include("$(pwd())/lib/poincareInitialData.jl")
include("$(pwd())/lib/ncmprojMKCONTROLGRID.jl")

function galerkinGUI()
	windowGalerkin = @Window("Galerkin Controls", 256, 256, false, true)
	setproperty!(windowGalerkin, :resizable, false)
	# setproperty!(windowGalerkin, :type_hint, Gtk.GdkWindowTypeHint.TOOLBAR)

	gridGalerkin = @Grid()
	setproperty!(gridGalerkin, :column_spacing, 5)
	setproperty!(gridGalerkin, :row_spacing, 5)
	push!(windowGalerkin, gridGalerkin)

	#=local initItSS = [
		("Trans. Iterations", Int, 0, 1e8, 1000),
		("SS Iterations", Int, 0, 1e8, 1000),
		("Trans. StepSize", Float64, .0, 1.0, .001),
		("SS StepSize", Float64, 0, 1.0, .001),
		("c₀", Float64, .0, 100.0, .1),
		("Periods", Int, 1, 128,1),
		("m", Int, 8, 4096, 1)
	]=#
	local initItSS = [
		("Trans. Iterations", Int, 0, 1e8, 1000),
		("Trans. StepSize", Float64, .0, 1.0, .001),
		("SS StepSize", Float64, 0, 1.0, .001),
		("Max. period", Int, 1, 128,1),
		("Intersections", Int, 1, 128,1),
		("m", Int, 8, 4096, 1),
		("c₀", Float64, .0, 100.0, .1)
	]
	dataItSS, gridItSS = mkControlGrid(initItSS, 1)
	gridGalerkin[1:2,1] = gridItSS

	buttonFindInitialValue = @Button("Find Initial Solution")
	setproperty!(buttonFindInitialValue, :expand, false)
	signal_connect(w -> @async(handlerFindInitialData(w, dataItSS)), buttonFindInitialValue, "clicked")
	gridGalerkin[1:2,2] = buttonFindInitialValue

	dataRS, gridRS = mkControlGrid([("Samples", Int, 8, 4096, 1)])
	buttonResample = @Button("Resample")
	setproperty!(buttonResample, :expand, false)
	signal_connect(w -> @async(handlerResample(w,dataRS)), buttonResample, "clicked")
	gridGalerkin[1,3], gridGalerkin[2,3] = gridRS, buttonResample

	showall(windowGalerkin)
end


function handlerResample(w, D)
	setproperty!(w, :sensitive, false)
	lock(lockProject)
	try
		ω,ℵ = project.activeSolution[end-1:end]
		f = toTrajInterp(project.activeSolution, 3)
		tmp = f(linspace(.0, 2pi, D["Samples"]))
		tmp = reduce(hcat, tmp)'
		tmp = rfft(tmp, [1])
		tmp = vec([real(tmp); imag(tmp[2:end,:])])
		tmp = [tmp; ω]

		Htmp(V) = H([V; ℵ])
		Jtmp(V) = J([V; ℵ])[:, 1:end-1]
		tmp = newton(Htmp, Jtmp, tmp, predCount(10) ∧ predEps(1e-10))
		project.activeSolution = [tmp; ℵ]
	finally
		unlock(lockProject)
		setproperty!(w, :sensitive, true)
	end
	return Void
end


function handlerFindInitialData(w, dataGUI)
	setproperty!(w, :sensitive, false)
	lock(lockProject)

	try
		#TODO function/macro bringIntoScope(D::Dict)
		# TIters, SSIters, TStepSize, SSStepSize, Periods, m, c₀ = map(x->dataGUI[x], ["Trans. Iterations", "SS Iterations", "Trans. StepSize", "SS StepSize", "Periods", "m", "c₀"])
		TIters, TStepSize, SSStepSize, maxCycles, nIntersections, m, c₀ = map(x->dataGUI[x], ["Trans. Iterations", "Trans. StepSize", "SS StepSize", "Max. period", "Intersections", "m", "c₀"])

		#=tmp = @fetch begin
			dataT, dataSS, P = findCycle((t,v)->f(t,[v;c₀]), .0, rand(3), TIters, TStepSize, SSIters, SSStepSize)
			cyc,ω = prepareCycle(dataSS, SSStepSize, P; fac=Periods)

			C = resample(rfft(cyc, [1]), m)
			C = [vec(vcat(real(C), imag(C[2:end,:]))); ω]

			Htmp(V) = H([V; c₀]) # R^N -> R^N system
			Jtmp(V) = J([V; c₀])[:, 1:end-1] # R^N -> R^(NxN) system
			return newton(Htmp, Jtmp, C, predCount(10) ∧ predEps(1e-10))
		end=#

		println("TIters: $TIters, TStepSize: $TStepSize, SSStepSize: $SSStepSize, maxCycles: $maxCycles, nIntersections: $nIntersections, m:$m, c₀:$c₀")

		tmp = @fetch begin
			cyc, P, T = findCyclePoincare((t,v)->f(t,[v;c₀]), rand(3),
				nIntersections=nIntersections, maxCycles=maxCycles, sampleSize=m,
				transientIterations=TIters, transientStepSize=TStepSize,
				steadyStateStepSize=SSStepSize)

			println("period $P for c=$(c₀) with T=$T on $(size(cyc,1)) samples")
			return

			ω = 1/T
			C = rfft(cyc, [1])
			C = [vec(vcat(real(C), imag(C[2:end,:]))); ω]

			Htmp(V) = H([V; c₀]) # R^N -> R^N system
			Jtmp(V) = J([V; c₀])[:, 1:end-1] # R^N -> R^(NxN) system
			return newton(Htmp, Jtmp, C, predCount(10) ∧ predEps(1e-10))
		end

		project.activeSolution = [tmp; c₀]
	finally
		unlock(lockProject)
		setproperty!(w, :sensitive, true)
	end
	return Void
end



function toTrajInterp(V, d)
	m = (length(V)-d-2)÷(2*d)
	tmp = reshape(V[1:end-2], 2m+1, d)
	return matsboUTIL.vectorize(x -> Float64[ matsboINTERPOLATE.interpolateTrigonometric(tmp[1,i], 2tmp[2:2+m-1,i], -2tmp[2+m:end,i])(x) / (2m+1) for i in 1:d ])
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


# takes real F-coefficients  V  , returns resampled version
function resample(V, m::Int)
	mapslices(V, [1]) do v
		f(x) = interpolateTrigonometric(real(v[1]), 2*real(v[2:end]), -2*imag(v[2:end]))(x) / (2*length(v)-1)
		rfft(f(linspace(.0,2pi,2*m+2)[1:end-1]))
	end
end
