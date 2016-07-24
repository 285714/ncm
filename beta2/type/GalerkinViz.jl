#TODO maintain projected state better
#		maybe forbid changing branch in middle?
#			=> consistency, can only add or remove point of projection!
#			callback from project abstraction
#TODO lock: lower opacity, hide solution dots

#TODO keep max/min, prevent relim
#TODO backspace from currentSolution

##### Global stuff
using PyCall
pygui(:tk) #prevent conflict with gtk
using PyPlot
D = PyCall.PyDict(matplotlib["rcParams"])

# try
# 	D["keymap.back"] =  deleteat!(
# 		D["keymap.back"],
# 		findfirst(D["keymap.back"], "backspace")
# 	)
# end

# remove all keybindings
for k in keys(D)
	contains(k, "keymap") && (D[k] = Any[])
end
ioff()

#####

# this visualization module expects functions specific to Galerkin (projection)

type GalerkinViz <: Visualization
	P::Project
	G::Galerkin
	fig
	idxLineBranch::Dict
	idxBranchLines::Dict
	activeSolutionMark
end

Base.show(V::GalerkinViz) = V.fig["show"]()

function GalerkinViz(P::Project, G::Galerkin)
	fig = figure()
	fig["canvas"]["set_window_title"]("Bifurcation Plot")

	V = GalerkinViz(P, G, fig, Dict{Any,Any}(), Dict{Any,Any}(), Void)

	fig["canvas"]["mpl_connect"]("pick_event", ev -> handlerPick(V,ev))
	fig["canvas"]["mpl_connect"]("key_press_event", ev -> handlerKey(V,ev))

	#NOTE events bound to project, not Viz obj! -> cannot change project!!!
	#connect to proj events
	observe(P, :pushBranch) do b
		pushBranch(V, b) end

	observe(P, :activeSolutionChanged) do
		s = V.P.activeSolution #TODO no direct access...
		s == Void && return
		typeof(s) == Solution && (s = s.data)
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

	# display the V.P
	map(b -> pushBranch(V,b), P)

	return V
end


function handlerPick(V::GalerkinViz, ev)
	b = get(V.idxLineBranch, ev[:artist], Void)
	b == Void && return Void
	try
		i = ev[:ind][1]+1
		setActiveSolution(V.P, b[i])
	end
	return Void
end

function handlerKey(V::GalerkinViz, ev)
	global tmp = ev #TODO asdf

	figure(V.fig[:number])

	#TODO implement dequeue style branch handling
	#TODO stop fooling with project internals
	#TODO make sure V.project is used, not global project
	#TODO locking...

	if ev[:key] == "backspace"
		!isa(V.P.activeSolution, Solution) && return Void
		s,b = V.P.activeSolution, V.P.activeSolution.parent
		if length(b) ≤ 2
			delBranch(V, b)
			deleteat!(V.P, findfirst(V.P, b)) #TODO fix
		elseif s==b[end]
			setActiveSolution(V.P, b[end-1])
			pop!(V.P.activeSolution.parent)
		elseif s==b[1]
			setActiveSolution(V.P, b[2])
			unshift!(V.P.activeSolution.parent)
		end
	elseif ev[:key] == "ctrl+delete" #TODO fix
		!isa(V.P.activeSolution, Solution) && return Void
		#delBranch(V, V.P.activeSolution.parent)
		b = V.P.activeSolution.parent
		delBranch(V, b)
		deleteat!(V.P, findfirst(V.P, b)) #TODO fix
	elseif ev[:key] == "r"
		ax = subplot(121)
		ax[:relim]()
		autoscale()
	elseif ev[:key] == "ctrl+r"
		clf()
		empty!(V.idxBranchLines)
		empty!(V.idxLineBranch)
		map(b -> pushBranch(V,b), V.P)
	elseif ev[:key] == "l"
		!isa(V.P.activeSolution, Solution) && return Void
		b = V.P.activeSolution.parent
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
		if typeof(V.P.activeSolution) == Solution
			for l in V.idxBranchLines[V.P.activeSolution.parent]
				l[:set_color](colMap[ev[:key]])
			end
		end
	end

	PyPlot.draw()
	return Void
end


function pushBranch(V::GalerkinViz, b::Branch)
	figure(V.fig[:number])
	subplot(121)

	x = map(last, b.solutions)
	y = reduce(map(transpose∘projection, b.solutions)) do A,a
		sA,sa = size(A,2), size(a,2)
		s = max(sA,sa)
		return [ A[:,(0:s-1)%sA+1]; a[1,(0:s-1)%sa+1] ]
	end

	lines = plot(x, y, picker=5, color="k", marker=".", markersize=3)
	PyPlot.draw()

	map(l -> (V.idxLineBranch[l]=b), lines)
	V.idxBranchLines[b] = lines

	return lines
end


function addToBranch(op, V::GalerkinViz, b::Branch, s::Solution)
	x = last(s) #TODO extend relay?
	Y = projection(s) #TODO broadcast to common size

	!haskey(V.idxBranchLines, b) && pushBranch(V,b)
	lines = V.idxBranchLines[b]

	for (l,y) in zip(lines, Y)
		l[:set_xdata](op(l[:get_xdata](), x))
		l[:set_ydata](op(l[:get_ydata](), y))
	end

	figure(V.fig[:number])
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

	figure(V.fig[:number])
	PyPlot.draw()
	return Void
end


# plotSolution(V::GalerkinViz, S::Solution) = plotSolution(B,S.data) # handled by convert? prob not..
function plotSolution(V::GalerkinViz, v::Vector{Float64})
	figure(V.fig[:number])

	subplot(121)
	try V.activeSolutionMark[:remove]() end
	x = v[end]
	y = projection(v)
	V.activeSolutionMark = scatter(fill(x, size(y)), y, facecolors="None", edgecolors="k", marker="o")

	subplot(122, projection="3d")
	t = linspace(0, 2pi, length(v)÷3 * 8)
	w = reduce(hcat, toTrajInterp(v, 3)(t))'
	plot(w[:,1], w[:,2], w[:,3], color="k")

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

	figure(V.fig[:number])
	PyPlot.draw()
end
