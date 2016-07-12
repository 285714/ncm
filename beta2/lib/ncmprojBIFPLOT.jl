# encapsulate plotting of bifurcation data
# provides a view on a project
# should be auto-updating for minimum interface

#TODO maintain projected state better
#		maybe forbid changing branch in middle?
#			=> consistency, can only add or remove point of projection!
#			callback from project abstraction
#TODO lock: lower opacity, hide solution dots

using PyCall
pygui(:tk) #prevent conflict with gtk
using PyPlot
D = PyCall.PyDict(matplotlib["rcParams"])
D["keymap.back"] =  deleteat!(
	D["keymap.back"],
	findfirst(D["keymap.back"], "backspace")
	)
ioff()

type BifPlot
	fig
	project::Project
	idxLineBranch::Dict
	idxBranchLines::Dict
	activeSolutionMark
end


function BifPlot(project::Project)
	fig = figure()
	fig["canvas"]["set_window_title"]("Bifurcation Plot")
	fig["show"]()

	#index for associating plotlines with solutions
	idxLineBranch = Dict{Any,Any}()
	idxBranchLines = Dict{Any,Any}()

	B = BifPlot(fig, project, idxLineBranch, idxBranchLines, Void)

	fig["canvas"]["mpl_connect"]("pick_event", ev -> handlerPick(B,ev))
	fig["canvas"]["mpl_connect"]("key_press_event", ev -> handlerKey(B,ev))

	#connect to project events
	observe(project, :pushBranch) do b
		pushBranch(B, b) end

	observe(project, :activeSolutionChanged) do
		s = project.activeSolution
		s == Void && return
		typeof(s) == Solution && (s = s.data)
		plotSolution(B, s)
	end

	observe(project, :pushSolution) do b,s
		addToBranch(push!, B, b, s) end

	observe(project, :unshiftSolution) do b,s
		addToBranch(unshift!, B, b, s) end

	observe(project, :popSolution) do b,s
		delFromBranch(pop!, B, b) end

	observe(project, :shiftSolution) do b,s
		delFromBranch(shift!, B, b) end

	observe(project, :delBranch) do b
		delBranch(B, b) end

	return B
end




function handlerPick(B::BifPlot, ev)
	b = get(B.idxLineBranch, ev[:artist], Void)
	b == Void && return Void
	lock(lockProject) #TODO move locking in this case to setActiveSolution
	try
		i = ev[:ind][1]+1
		setActiveSolution(project, b[i])
	finally
		unlock(lockProject)
	end
	return Void
end

function handlerKey(B::BifPlot, ev)
	global tmp = ev #TODO asdf

	figure(B.fig[:number])

	#TODO implement dequeue style branch handling
	#TODO stop fooling with project internals
	#TODO make sure B.project is used, not global project
	#TODO locking...

	if ev[:key] == "delete"
		!isa(B.project.activeSolution, Solution) && return Void
		s,b = B.project.activeSolution, B.project.activeSolution.parent
		length(b) > 1 && s==b[1] && setActiveSolution(B.project, b[2])
		shift!(B.project.activeSolution.parent)
	elseif ev[:key] == "backspace"
		!isa(B.project.activeSolution, Solution) && return Void
		s,b = B.project.activeSolution, B.project.activeSolution.parent
		length(b) > 1 && s==b[end] && setActiveSolution(B.project, b[end-1])
		pop!(B.project.activeSolution.parent)
	elseif ev[:key] == "ctrl+delete" #TODO fix
		!isa(B.project.activeSolution, Solution) && return Void
		#delBranch(B, B.project.activeSolution.parent)
		b = B.project.activeSolution.parent
		deleteat!(B.project, findfirst(B.project, b)) #TODO fix
		delBranch(B, b)
	elseif ev[:key] == "r"
		ax = subplot(121)
		ax[:relim]()
		autoscale()
	elseif ev[:key] == "ctrl+r"
		clf()
		empty!(B.idxBranchLines)
		empty!(B.idxLineBranch)
		map(b -> pushBranch(B,b), B.project)
	elseif ev[:key] == "ctrl+l"
		!isa(B.project.activeSolution, Solution) && return Void
		b = B.project.activeSolution.parent
		for l in B.idxBranchLines[b]
			if l[:get_picker]() != 5
				l[:set_picker](5)
				l[:set_marker](".")
			else
				l[:set_picker](0)
				l[:set_marker]("None")
			end
		end
	elseif ev[:key] == "ctrl+L"
		for l in keys(B.idxLineBranch)
			l[:set_picker](5)
			l[:set_marker](".")
		end
	elseif ev[:key] == "ctrl+alt+L"
		for l in keys(B.idxLineBranch)
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
		if typeof(project.activeSolution) == Solution
			for l in B.idxBranchLines[project.activeSolution.parent]
				l[:set_color](colMap[ev[:key]])
			end
		end
	end

	PyPlot.draw()
	return Void
end


function pushBranch(B::BifPlot, b::Branch)
	figure(B.fig[:number])
	subplot(121)

	x = map(last, b.solutions)
	y = reduce(map(transpose∘projection, b.solutions)) do A,a
		sA,sa = size(A,2), size(a,2)
		s = max(sA,sa)
		return [ A[:,(0:s-1)%sA+1]; a[1,(0:s-1)%sa+1] ]
	end

	lines = plot(x, y, picker=5, color="k", marker=".", markersize=3)
	PyPlot.draw()

	map(l -> (B.idxLineBranch[l]=b), lines)
	B.idxBranchLines[b] = lines

	return lines
end


function addToBranch(op, B::BifPlot, b::Branch, s::Solution)
	x = last(s) #TODO extend relay?
	Y = projection(s) #TODO broadcast to common size

	lines = B.idxBranchLines[b]
	for (l,y) in zip(lines, Y)
		l[:set_xdata](op(l[:get_xdata](), x))
		l[:set_ydata](op(l[:get_ydata](), y))
	end

	figure(B.fig[:number])
	PyPlot.draw()
	return Void
end

function delFromBranch(op, B::BifPlot, b::Branch)
	for l in B.idxBranchLines[b]
		x,y = l[:get_xdata](), l[:get_ydata]()
		op(x); op(y)
		l[:set_xdata](x); l[:set_ydata](y)
	end

	figure(B.fig[:number])
	PyPlot.draw()
	return Void
end


# plotSolution(B::BifPlot, S::Solution) = plotSolution(B,S.data) # handled by convert? prob not..
function plotSolution(B::BifPlot, V::Vector{Float64})
	figure(B.fig[:number])

	subplot(121)
	try B.activeSolutionMark[:remove]() end
	x = V[end]
	y = projection(V)
	B.activeSolutionMark = scatter(fill(x, size(y)), y, facecolors="None", edgecolors="k", marker="o")

	subplot(122, projection="3d")
	t = linspace(0, 2pi, length(V)÷3 * 8)
	v = reduce(hcat, toTrajInterp(V, 3)(t))'
	p2 = plot(v[:,1], v[:,2], v[:,3], color="k")

	PyPlot.draw()
	return Void
end


function delBranch(B::BifPlot, b::Branch)
	lines = B.idxBranchLines[b]
	delete!(B.idxBranchLines, b)
	map(lines) do l
		l[:remove]()
		delete!(B.idxLineBranch, l)
	end

	figure(B.fig[:number])
	PyPlot.draw()
end
