using PyCall
pygui(:tk) #prevent conflict with gtk
using PyPlot

using Gtk.ShortNames
using matsboINTERPOLATE, matsboNWTN, matsboPRED, matsboUTIL

include("../lib/ncmprojFINDINITIALDATA.jl")
include("../lib/ncmprojMKCONTROLGRID.jl")

function systemExec()
	roesslerGUI()

	global bifurcationFig = figure()
	bifurcationFig["canvas"]["set_window_title"]("Bifurcation Plot")
	bifurcationFig["canvas"]["mpl_connect"]("pick_event", handlerBifurcationPick)

	#index for associating plotlines with solutions
	global idxLinesToBranch = Dict{PyObject, Int}()

	global solutionFig = figure()
	solutionFig["canvas"]["set_window_title"]("Solution Plot")
end


function handlerBifurcationPick(ev)
	i = idxLinesToBranch[ev["artist"]]
	j = ev["ind"][1]+1
	plotBranch(project.branches[i])
	plotSolution(project.branches[i].solutions[j])
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
		("SS Iterations", Int, 0, 1e8, 1000),
		("Trans. StepSize", Float64, .0, 1.0, .001),
		("SS StepSize", Float64, 0, 1.0, .001),
		("c₀", Float64, .0, 100.0, .1),
		("Periods", Int, 1, 128,1),
		("m", Int, 0, 4096, 1)
	]
	dataItSS, gridItSS = mkControlGrid(initItSS, 1)
	gridRoessler[1,1] = gridItSS

	buttonFindInitialValue = @Button("Find Initial Value")
	signal_connect(w -> @async(handlerFindInitialData(dataItSS)), buttonFindInitialValue, "clicked")
	gridRoessler[1,2] = buttonFindInitialValue

	showall(windowRoessler)
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
	plotSolution(C)
end


function plotSolution(V)
	f = toTrajInterp(V, 3)
	t = linspace(.0,2pi,1024)
	T = reduce(hcat, f(t))'

	global solutionFig
	figure(solutionFig[:number])
	clf()
	gca(projection="3d")
	plot(T[:,1], T[:,2], T[:,3])

	plotBifurcationSingle(V)
end


function projection(V)
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


function plotBifurcation(project::Project)
	global bifurcationFig
	figure(bifurcationFig[:number])
	clf()

	#TODO parallel, changing length of projection
	for i in 1:length(project.branches)
		branch = project.branches[i]

		if true #get(branch.D, "changed", true)
			try branch.D["plot"]["remove"]() end
			x = map(last, branch.solutions)
			y = reduce(hcat, map(projection, branch.solutions))'
			lines = plot(x,y,picker=5)
			map(l -> (idxLinesToBranch[l]=i), lines)
			branch.D["plot"] = lines
			branch.D["changed"] = false
		end
	end
end


function plotBranch(branch::Branch)
	global bifurcationFig, branchPoints
	figure(bifurcationFig[:number])

	try for b in branchPoints try b["remove"]() end end end
	branchPoints = Any[]

	x = map(last, branch.solutions)
	y = reduce(vcat, map(transpose∘projection, branch.solutions))
	scatter(broadcast((x,y)->x, x, y), y, color="k")

	return Void
end


function plotBifurcationSingle(V)
	global bifurcationFig, singleSolutionMark
	figure(bifurcationFig[:number])

	try singleSolutionMark["remove"]() end

	y = projection(V)
	x = fill(V[end], size(y))
	singleSolutionMark = scatter(x,y,color="r")

	return Void
end



function Hroessler(V::Vector{Float64})
	local m = length(V-5)÷6
	local X₀,Xᵣ,Xᵢ,Y₀,Yᵣ,Yᵢ,Z₀,Zᵣ,Zᵢ,ω,ℵ

	ω,ℵ = V[end-1:end]
	V = reshape(V[1:end-2],2m+1,3)
	X₀,Y₀,Z₀ = V[1,1], V[1,2], V[1,3]
	Xᵣ,Yᵣ,Zᵣ = V[2:2+m-1,1], V[2:2+m-1,2], V[2:2+m-1,3]
	Xᵢ,Yᵢ,Zᵢ = V[2+m:end,1], V[2+m:end,2], V[2+m:end,3]

	local S₀,Sᵣ,Sᵢ
	local tmp = rfft( irfft(complex([X₀;Xᵣ], [.0;Xᵢ]), 2m+1) .* irfft(complex([Z₀;Zᵣ], [.0;Zᵢ]), 2m+1) )
	S₀,Sᵣ,Sᵢ = real(tmp[1]), real(tmp[2:end]), imag(tmp[2:end])

	local D = ω*(1:m)

	rtn = [
		#A
		-Y₀-Z₀
		-Yᵣ-Zᵣ+D.*Xᵢ
		-Yᵢ-Zᵢ-D.*Xᵣ
		#B
		X₀+a*Y₀
		Xᵣ+a*Yᵣ+D.*Yᵢ
		Xᵢ+a*Yᵢ-D.*Yᵣ
		#C
		-ℵ*Z₀+S₀+(2m+1)*b
		-ℵ*Zᵣ+Sᵣ+D.*Zᵢ
		-ℵ*Zᵢ+Sᵢ-D.*Zᵣ

		X₀+sum(Xᵣ)
	]

  return rtn
end



function Jroessler(V)
	local m = length(V-5)÷6
	local X₀,Xᵣ,Xᵢ,Y₀,Yᵣ,Yᵢ,Z₀,Zᵣ,Zᵢ,ω,ℵ

	ω,ℵ = V[end-1:end]
	V = reshape(V[1:end-2],2m+1,3)
	X₀,Y₀,Z₀ = V[1,1], V[1,2], V[1,3]
	Xᵣ,Yᵣ,Zᵣ = V[2:2+m-1,1], V[2:2+m-1,2], V[2:2+m-1,3]
	Xᵢ,Yᵢ,Zᵢ = V[2+m:end,1], V[2+m:end,2], V[2+m:end,3]

	local D = ω*(1:m)

	# Jacobian of circular convolution of coefficients of real functions in appropriate format.
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

	local AbyX, AbyY, AbyZ, Abyω, Abyℵ,
		BbyX, BbyY, BbyZ, Bbyω, Bbyℵ,
		CbyX, CbyY, CbyZ, Cbyω, Cbyℵ,
		LbyX, LbyY, LbyZ, Lbyω, Lbyℵ

	AbyX = [
		.0					zeros(1,m)			zeros(1,m)
		zeros(m)			zeros(m,m)			diagm(D)
		zeros(m)			diagm(-D)			zeros(m,m)
	]

	AbyY = AbyZ = -eye(2m+1)

	BbyX = eye(2m+1)

	BbyY = [
		a					zeros(1,m)			zeros(1,m)
		zeros(m)			a*eye(m)			diagm(D)
		zeros(m)			diagm(-D)			a*eye(m)
	]

	BbyZ = zeros(2m+1,2m+1)

	CbyX = rCCD(Z₀,Zᵣ,Zᵢ)

	CbyY = zeros(2m+1, 2m+1)

	CbyZ = rCCD(X₀,Xᵣ,Xᵢ) + [
		-ℵ					zeros(1,m)			zeros(1,m)
		zeros(m)			-ℵ*eye(m,m)			diagm(D)
		zeros(m)			diagm(-D)			-ℵ*eye(m,m)
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

	return [
		AbyX AbyY AbyZ Abyω Abyℵ
		BbyX BbyY BbyZ Bbyω Bbyℵ
		CbyX CbyY CbyZ Cbyω Cbyℵ
		LbyX LbyY LbyZ Lbyω Lbyℵ
	]
end


H,J = Hroessler, Jroessler




function toTrajInterp(V, d)
	m = (length(V)-d-2)÷(2*d)
	tmp = reshape(V[1:end-2], 2m+1, d)
	return matsboUTIL.vectorize(x -> [ matsboINTERPOLATE.interpolateTrigonometric(tmp[1,i], 2tmp[2:2+m-1,i], -2tmp[2+m:end,i])(x) / (2m+1) for i in 1:d ])
end
