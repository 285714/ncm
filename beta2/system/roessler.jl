using PyPlot
include("../lib/ncmprojFINDINITIALDATA.jl")
include("../lib/ncmprojMKCONTROLGRID.jl")

function systemExec()
	roesslerGUI()

	global bifurcationFig = figure(figsize=(9,5))
	bifurcationFig["canvas"]["set_window_title"]("Bifurcation Plot")

	global solutionFig = figure(figsize=(9,5))
	solutionFig["canvas"]["set_window_title"]("Solution Plot")
end



a,b = .1,.1
roessler(v) = [-v[2] - v[3]; v[1] + a*v[2]; b + v[3]*(v[1]-v[4])]
roessler(t,v) = roessler(v)

function roesslerGUI()
	windowRoessler = @Window("Roessler Controls", 128, 256, false, true)
	# setproperty!(windowGalerkin, :type_hint, Gtk.GdkWindowTypeHint.TOOLBAR)

	gridRoessler = @Grid()
	setproperty!(gridRoessler, :column_spacing, 5)
	setproperty!(gridRoessler, :row_spacing, 5)
	push!(windowRoessler, gridRoessler)

	local dataItSS, gridItSS
	local initItSS = [
		("Trans. Iterations", Int, 0, 1e8, 1000),
		("Trans. StepSize", Float64, :h, .0, 1.0, .001),
		("m", Int, 0, 4096, 1),
		("SS Iterations", Int, 0, 1e8, 1000),
		("SS StepSize", Float64, :h, 0, 1.0, .001),
		("Periods", Int, 1, 128,1)
	]
	dataItSS, gridItSS = mkControlGrid(initItSS, 2)
	gridRoessler[1,1] = gridItSS

	buttonFindInitialValue = @Button("Find Initial Value")
	signal_connect(w -> handlerFindInitialData(), buttonFindInitialValue, "clicked")
	gridRoessler[1,2] = buttonFindInitialValue

	showall(windowRoessler)
end



function handlerFindInitialData()
	dataT, dataSS, P = findCycle((t,v)->f(t,[v;7.0]), .0, rand(d), TIters, TStepSize, SSIters, SSStepSize)
	cyc,ω = prepareCycle(dataSS, SSStepSize, P; fac=Periods)

	C = mapslices(cyc, [1]) do v
		tmp = rfft(v)
		m = length(tmp)-1
		f = x -> interpolateTrigonometric(real(tmp[1]), 2*real(tmp[2:end]), -2*imag(tmp[2:end]))(x) / (2m+1)
		tmp = map(f, linspace(.0,2pi,2*m+2)[1:end-1])
		tmp = rfft(tmp)
		return [ real(tmp); imag(tmp[2:end]) ]
	end

	C = [reduce(vcat, C); ω; 7.0]
end

function plotSolution(C)
	global solutionFig
	plot()
end
