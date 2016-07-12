#TODO encapsulate... project, dimensions, data

using Gtk.ShortNames
using mbInterpolate, mbNewton, mbPred, mbUtil
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

	dataGal = Dict{AbstractString, Any}()
	gridGalerkin[1:2,1] = mkControlGrid(dataGal, [
		("Trans. Iterations", Int, 2000, 1000:1000:1e8),
		("Trans. StepSize", Float64, .1, 1e-3:.001:1.0),
		("SS StepSize", Float64, .1, 1e-3:.001:1.0),
		("Max. period", Int, 30, 1:128),
		("Intersections", Int, 120, 1:128),
		("m", Int, 64, 8:4096),
		("c₀", Float64, .0, -500:.01:500)
	])

	buttonFindInitialValue = @Button("Find Initial Solution")
	setproperty!(buttonFindInitialValue, :expand, false)
	signal_connect(w -> @async(handlerFindInitialData(w, dataGal)), buttonFindInitialValue, "clicked")
	gridGalerkin[1:2,2] = buttonFindInitialValue

	gridGalerkin[1:2,3] = mkControlGrid(dataGal, [
		("Samples", Int, 64, 8:4096),
		("Crop", Float64, 1.0, .0:.1:8.0)
	])
	buttonResample = @Button("Process")
	setproperty!(buttonResample, :expand, false)
	signal_connect(w -> @async(handlerResample(w, dataGal)), buttonResample, "clicked")
	gridGalerkin[1:2,4] = buttonResample

	showall(windowGalerkin)
end


function handlerResample(w, D)
	setproperty!(w, :sensitive, false)
	lock(project)
	try
		ω,ℵ = project.activeSolution[end-1:end]
		m = (length(project.activeSolution)-5)÷6
		V = reshape(project.activeSolution[1:end-2], 2m+1, 3)
		V = complex(V[1:m+1,:], [zeros(1, size(V,2)); V[m+2:end,:]])

		C = resample(V, D["Samples"], D["Crop"])
		C = [vec(vcat(real(C), imag(C[2:end,:]))); ω/D["Crop"]]

		Htmp(V) = H([V; ℵ])
		Jtmp(V) = J([V; ℵ])[:, 1:end-1]
		tmp = newton(Htmp, Jtmp, C, predCount(10) ∧ predEps(1e-10))

		setActiveSolution(project, [tmp; ℵ])
	finally
		unlock(project)
		setproperty!(w, :sensitive, true)
	end
	return Void
end


function handlerFindInitialData(w, dataGal)
	setproperty!(w, :sensitive, false)
	lock(project)

	try
		#TODO function/macro bringIntoScope(D::Dict)
		# TIters, SSIters, TStepSize, SSStepSize, Periods, m, c₀ = map(x->dataGal[x], ["Trans. Iterations", "SS Iterations", "Trans. StepSize", "SS StepSize", "Periods", "m", "c₀"])
		TIters, TStepSize, SSStepSize, maxCycles, nIntersections, m, c₀ = map(x->dataGal[x], ["Trans. Iterations", "Trans. StepSize", "SS StepSize", "Max. period", "Intersections", "m", "c₀"])

		tmp = @fetch begin
			# dataT, dataSS, P = findCycle((t,v)->f(t,[v;c₀]), .0, rand(3), TIters, TStepSize, SSIters, SSStepSize)
			# cyc,ω = prepareCycle(dataSS, SSStepSize, P; fac=Periods)
			cyc, P, ω = findCyclePoincare((t,v)->f(t,[v;c₀]), rand(3),
				nIntersections=nIntersections, maxCycles=maxCycles, sampleSize=m,
				transientIterations=TIters, transientStepSize=TStepSize,
				steadyStateStepSize=SSStepSize)

			C = resample(rfft(cyc, [1]), m)
			C = [vec(vcat(real(C), imag(C[2:end,:]))); ω]

			Htmp(V) = H([V; c₀]) # R^N -> R^N system
			Jtmp(V) = J([V; c₀])[:, 1:end-1] # R^N -> R^(NxN) system
			return newton(Htmp, Jtmp, C, predCount(10) ∧ predEps(1e-10))
		end

		setActiveSolution(project, [tmp; c₀])
	finally
		unlock(project)
		setproperty!(w, :sensitive, true)
	end
	return Void
end



function toTrajInterp(V, d)
	m = (length(V)-d-2)÷(2*d)
	tmp = reshape(V[1:end-2], 2m+1, d)
	return mbUtil.vectorize(x -> Float64[ mbInterpolate.interpolateTrigonometric(tmp[1,i], 2tmp[2:2+m-1,i], -2tmp[2+m:end,i])(x) / (2m+1) for i in 1:d ])
end


function projection(V)
	f = toTrajInterp(V,3)
	rtn = Float64[]

	dt = 2pi/1025 #TODO ...
	for t in .0:dt:2pi
		if f(t)[1] ≤ .0 ≤ f(t+dt)[1]
			x = mbUtil.bisection(x->f(x)[1], t, t+dt, ϵ=1e-4) #TODO precision
			push!(rtn, norm(f(x)))
		end
	end

	return rtn
end


# takes real F-coefficients  V  , returns resampled version
function resample(V, m::Int, scale=1.0)
	mapslices(V, [1]) do v
		f(x) = mbInterpolate.interpolateTrigonometric(real(v[1]), 2*real(v[2:end]), -2*imag(v[2:end]))(x) / (2*length(v)-1)
		rfft(f(linspace(.0,2pi*scale,2*m+2)[1:end-1]))
	end
end
