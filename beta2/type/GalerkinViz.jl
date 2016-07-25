#TODO keep max/min, prevent relim

##### Global stuff
using PyCall
pygui(:tk) #prevent conflict with gtk
using PyPlot
# remove all keybindings
D = PyCall.PyDict(matplotlib["rcParams"])
for k in keys(D); contains(k, "keymap") && (D[k] = Any[]) end
ioff()

#####

# this visualization module expects functions specific to Galerkin (projection)

type GalerkinViz <: Visualization
	parent::Session
	figBif
	figSol
	idxLineBranch::Dict
	idxBranchLines::Dict
	activeSolutionMark
end

function Base.show(V::GalerkinViz)
	figure(V.figBif[:number])
	clf()
	empty!(V.idxBranchLines)
	empty!(V.idxLineBranch)
	map(b -> pushBranch(V,b), V.parent.P)
	V.figBif["show"]()
	V.figSol["show"]()
end

function GalerkinViz(parent::Session)
	figBif = figure(figsize=(8,6), dpi=80, facecolor="w")
	figBif["canvas"]["set_window_title"]("Bifurcation Plot")
	subplot(111)
	subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

	figSol = figure(figsize=(4,3), dpi=80, facecolor="w")
	figSol["canvas"]["set_window_title"]("Solution Plot")
	subplot(111, projection="3d")
	subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
	figSol[:canvas][:toolbar][:pack_forget]()

	V = GalerkinViz(parent, figBif, figSol, Dict{Any,Any}(), Dict{Any,Any}(), Void)

	figBif["canvas"]["mpl_connect"]("pick_event", ev -> handlerPick(V,ev))
	figBif["canvas"]["mpl_connect"]("key_press_event", ev -> @schedule handlerKey(V,ev))

	#connect to proj events
	P = parent.P

	observe(P, :pushBranch) do b
		pushBranch(V, b) end

	observe(P, :activeSolutionChanged) do
		s = V.parent.P.activeSolution #TODO no direct access...
		!isa(s, Solution) && !isa(s, Vector{Float64}) && return
		isa(s, Solution) && (s = s.data)
		plotSolution(V, s)
	end

	observe(P, :pushSolution) do b,s
		addToBranch(push!, V, b, s) end

	observe(P, :unshiftSolution) do b,s
		addToBranch(unshift!, V, b, s) end

	observe(P, :popSolution) do b,s
		delFromBranch(pop!, V, b) end

	observe(P, :shiftSolution) do b,s
		delFromBranch(shift!, V, b) end

	observe(P, :delBranch) do b
		delBranch(V, b) end

	return V
end


function handlerPick(V::GalerkinViz, ev)
	b = get(V.idxLineBranch, ev[:artist], nothing)
	b == nothing && return

	try
		i = ev[:ind][1]+1
		setActiveSolution(V.parent.P, b[i])
	catch e
		println(e)
	end
	return Void
end

function handlerKey(V::GalerkinViz, ev)
	# println(ev[:key]) #debug

	figure(V.figBif[:number])

	#TODO implement dequeue style branch handling
	#TODO stop fooling with project internals
	#TODO make sure V.project is used, not global project
	#TODO locking...

	if ev[:key] == "backspace"
		!isa(V.parent.P.activeSolution, Solution) && return Void
		s,b = V.parent.P.activeSolution, V.parent.P.activeSolution.parent
		if length(b) ≤ 2
			delBranch(V, b)
			deleteat!(V.parent.P, findfirst(V.parent.P, b)) #TODO fix
		elseif s==b[end]
			setActiveSolution(V.parent.P, b[end-1])
			pop!(V.parent.P.activeSolution.parent)
		elseif s==b[1]
			setActiveSolution(V.parent.P, b[2])
			shift!(V.parent.P.activeSolution.parent)
		end
	elseif ev[:key] == "ctrl+delete" #TODO fix
		!isa(V.parent.P.activeSolution, Solution) && return Void
		b = V.parent.P.activeSolution.parent
		delBranch(V, b)
		deleteat!(V.parent.P, findfirst(V.parent.P, b)) #TODO fix
	elseif ev[:key] == "r"
		figure(V.figBif[:number])
		gca()[:relim]()
		autoscale()
	elseif ev[:key] == "ctrl+r"
		figure(V.figBif[:number])
		clf()
		empty!(V.idxBranchLines)
		empty!(V.idxLineBranch)
		map(b -> pushBranch(V,b), V.parent.P)
	elseif ev[:key] == "l"
		!isa(V.parent.P.activeSolution, Solution) && return Void
		b = V.parent.P.activeSolution.parent
		for l in V.idxBranchLines[b]
			if l[:get_picker]() != 5
				l[:set_picker](5)
				l[:set_marker](".")
			else
				l[:set_picker](0)
				l[:set_marker]("None")
			end
		end
	elseif ev[:key] == "ctrl+l"
		for l in keys(V.idxLineBranch)
			l[:set_picker](5)
			l[:set_marker](".")
		end
	elseif ev[:key] == "ctrl+L"
		for l in keys(V.idxLineBranch)
			l[:set_picker](0)
			l[:set_marker]("None")
		end
	elseif ev[:key] in [ "$(i)" for i in 0:9 ]
		const colMap = Dict(
			"0" => "#000000", "1" => "#FF0000",
			"2" => "#00FF00", "3" => "#0000FF",
			"4" => "#00FFFF", "5" => "#FF00FF",
			"6" => "#FFFF00", "7" => "#800000",
			"8" => "#008000", "9" => "#000080"
			)
		if typeof(V.parent.P.activeSolution) == Solution
			for l in V.idxBranchLines[V.parent.P.activeSolution.parent]
				l[:set_color](colMap[ev[:key]])
			end
		end
	elseif ev[:key] == "b"
		!isa(V.parent.P.activeSolution, Solution) && return
		b = V.parent.P.activeSolution.parent
		plotBranch(V, b)
	elseif ev[:key] == " "
		step(V.parent.cont)
	elseif ev[:key] == "e"
		!isa(V.parent.P.activeSolution, Solution) && return
		s = V.parent.P.activeSolution
		b = s.parent
		setActiveSolution(V.parent.P, s==b[1] ? b[end] : b[1])
	elseif ev[:key] in ["up","down"]
		P = V.parent.P
		s = P.activeSolution
		!isa(s, Solution) && return
		b = s.parent
		i = findfirst(P,b)-1
		b′ = P[ mod( ev[:key] == "up" ? i-1 : i+1, length(P)) + 1 ]
		setActiveSolution(P, b′[1])
	elseif ev[:key] in ["left","right"]
		P = V.parent.P
		s = P.activeSolution
		!isa(s, Solution) && return
		b = s.parent
		i = findfirst(b,s)-1
		s′ = b[ mod( ev[:key] == "left" ? i-1 : i+1 , length(b)) + 1 ]
		setActiveSolution(P, s′)
	end

	PyPlot.draw()
	return Void
end


function pushBranch(V::GalerkinViz, b::Branch)
	#TODO catch empty
	x = map(last, b)

	function reducer(A,a)
		sA,sa = size(A,2), size(a,2)
		s = max(sA,sa)
		return [
			A[:,(0:s-1)%sA+1]
			a[:,(0:s-1)%sa+1]
			]
	end

	mapper(x) = map(norm, projection(x))'

	y = mapreduce(mapper, reducer, b)

	figure(V.figBif[:number])
	hold(true)
	lines = plot(x, y, picker=5, color="k", marker=".", markersize=3)
	PyPlot.draw()

	map(l -> (V.idxLineBranch[l]=b), lines)
	V.idxBranchLines[b] = lines

	return lines
end


function addToBranch(op, V::GalerkinViz, b::Branch, s::Solution)
	x = last(s) #TODO extend relay?
	Y = map(norm, projection(s)) #TODO broadcast to common size

	!haskey(V.idxBranchLines, b) && error("branch missing")
	# !haskey(V.idxBranchLines, b) && pushBranch(V,b)
	lines = V.idxBranchLines[b]

	for (l,y) in zip(lines, Y)
		l[:set_xdata](op(l[:get_xdata](), x))
		l[:set_ydata](op(l[:get_ydata](), y))
	end

	figure(V.figBif[:number])
	PyPlot.draw()
	return Void
end

#TODO handle last solution removed?
function delFromBranch(op, V::GalerkinViz, b::Branch)
	for l in V.idxBranchLines[b]
		x,y = l[:get_xdata](), l[:get_ydata]()
		op(x); op(y)
		l[:set_xdata](x); l[:set_ydata](y)
	end

	figure(V.figBif[:number])
	PyPlot.draw()
	return Void
end


# plotSolution(V::GalerkinViz, S::Solution) = plotSolution(B,S.data) # handled by convert? prob not..
function plotSolution(V::GalerkinViz, v::Vector{Float64})
	intersections = projection(v)

	figure(V.figBif[:number])
	hold(true)
	try V.activeSolutionMark[:remove]() end
	x = v[end]
	y = map(norm, intersections)
	V.activeSolutionMark = scatter(fill(x, size(y)), y; facecolors="None", edgecolors="k", marker="o")
	PyPlot.draw()

	#TODO mark proj points

	figure(V.figSol[:number])
	cla()
	gca(projection="3d")
	hold(false)
	t = linspace(0, 2pi, length(v)÷3 * 8)
	w = reduce(hcat, toTrajInterp(v, 3)(t))'
	plot(w[:,1], w[:,2], w[:,3], color="k") #+ linspace(0, maximum(w[:,3]), length(w[:,3]))
	u = reduce(hcat, intersections)
	hold(true)
	scatter3D(u[1,:], u[2,:], u[3,:]; facecolors="None", edgecolors="k", marker="o")
	PyPlot.draw()

	return Void
end


function delBranch(V::GalerkinViz, b::Branch)
	!haskey(V.idxBranchLines, b) && return
	lines = V.idxBranchLines[b]
	delete!(V.idxBranchLines, b)
	map(lines) do l
		l[:remove]()
		delete!(V.idxLineBranch, l)
	end

	figure(V.figSol[:number])
	PyPlot.draw()
end


# plot a whole branch in
function plotBranch(V::GalerkinViz, B::Branch)
	figure(V.figSol[:number])
	clf()
	gca(projection="3d")
	hold(true)
	for v in B
		t = linspace(0, 2pi, length(v)÷3 * 8)
		w = reduce(hcat, toTrajInterp(v, 3)(t))'
		plot(w[:,1], w[:,2], w[:,3], color="k", alpha=.1)
	end
	PyPlot.draw()
	return
end
