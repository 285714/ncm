# encapsulate plotting of bifurcation data
# provides a view on a project
# should be auto-updating for minimum interface

#TODO maintain projected state better
#		maybe forbid changing branch in middle?
#			=> consistency, can only add or remove point of projection!
#			callback from project abstraction

using PyCall
pygui(:tk) #prevent conflict with gtk
using PyPlot
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

	function handlerPick(ev)
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

	function handlerKey(ev)
		global tmp = ev

		figure(B.fig[:number])

		#TODO implement dequeue style branch handling

		if ev[:key] == "delete"
			project.activeSolution == Void && return Void
			i,j = findSolution(project, project.activeSolution) #TODO fix this
			if i ≠ 0 && j ≠ 0
				deleteat!(project.branches[i].solutions, j)
				isempty(project.branches[i].solutions) && deleteat!(project.branches, i)
			end
		elseif ev[:key] == "ctrl+delete"
			S = project.activeSolution
			(S == Void || typeof(S) != Solution) && return Void
			b = S.parent
			lines = B.idxBranchLines[b]
			delete!(B.idxBranchLines, b)
			for l in lines
				delete!(B.idxLineBranch,l)
				l[:remove]()
			end
			deleteat!(b.parent.branches, findfirst(b.parent.branches, S))
		elseif ev[:key] == "r"
			ax = subplot(121)
			ax[:relim]()
			autoscale()
		# elseif ev[:key] == "ctrl+r"
		# 	empty!(B.cache)
		# 	empty!(B.idx)
		# 	clf()
		elseif ev[:key] in [ "$(i)" for i in 0:9 ]
			const colMap = Dict(
				"0" => "#000000", "1" => "#FF0000",
				"2" => "#00FF00", "3" => "#0000FF",
				"4" => "#00FFFF", "5" => "#FF00FF",
				"6" => "#FFFF00", "7" => "#800000",
				"8" => "#008000", "9" => "#000080"
				)
			try map(x->x[:set_color](colMap[ev[:key]]), B.currPlot) end
		end

		PyPlot.draw()

		return Void
	end

	fig["canvas"]["mpl_connect"]("pick_event", handlerPick)
	fig["canvas"]["mpl_connect"]("key_press_event", handlerKey)

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
		addToBranch(push!, B::BifPlot, b::Branch, s::Solution) end

	observe(project, :unshiftSolution) do b,s
		addToBranch(unshift!, B::BifPlot, b::Branch, s::Solution) end

	observe(project, :popSolution) do b,s
		# update plot
	end

	observe(project, :shiftSolution) do b,s
		# update plot
	end

	return B
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

function delFromBranch(op, B::BifPlot, b::Branch, s::Solution)
	error("asdfasdfasdf")
end


# plotSolution(B::BifPlot, S::Solution) = plotSolution(B,S.data) # handled by convert? prob not..
function plotSolution(B::BifPlot, V::Vector{Float64})
	figure(B.fig[:number])

	subplot(121)
	try B.activeSolutionMark[:remove]() end
	x = V[end]
	y = projection(V)
	B.activeSolutionMark = scatter(fill(x, size(y)), y, color="r")

	subplot(122, projection="3d")
	t = linspace(0, 2pi, length(V)÷3 * 8)
	v = reduce(hcat, toTrajInterp(V, 3)(t))'
	p2 = plot(v[:,1], v[:,2], v[:,3], color="k")

	PyPlot.draw()

	return Void
end
