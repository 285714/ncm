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
	idx::Dict
	cache::Dict
	solutionMark
end


function BifPlot(project::Project)
	fig = figure()
	fig["canvas"]["set_window_title"]("Bifurcation Plot")
	fig["show"]()

	#index for associating plotlines with solutions
	idx = Dict{Any,Any}()
	cache = Dict{Any,Any}()

	B = BifPlot(fig, project, idx, cache, Void)

	function handlerPick(ev)
		branch = get(B.idx, ev["artist"], Void)
		branch == Void && return Void
		i = ev["ind"][1]+1
		# lock(lockProject)
		project.activeSolution = branch.solutions[i]
		# unlock(lockProject)
		# B.cache[project.activeSolution] = plotSolution(B, project.activeSolution)
		return Void
	end

	function handlerKey(ev)
		global tmp = ev
		if ev[:key] == "delete"
			project.activeSolution == Void && return Void
			i,j = findSolution(project, project.activeSolution)
			if i ≠ 0 && j ≠ 0
				deleteat!(project.branches[i].solutions, j)
				isempty(project.branches[i].solutions) && deleteat!(project.branches, i)
			end
		elseif ev[:key] == "ctrl+delete"
			project.activeSolution == Void && return Void
			i,j = findSolution(project, project.activeSolution)
			i ≠ 0 && deleteat!(project.branches, i)
		elseif ev[:key] == "r"
			figure(B.fig[:number])
			ax = subplot(121)
			ax[:relim]()
			autoscale()
		elseif ev[:key] == "ctrl+r"
			empty!(B.cache)
			empty!(B.idx)
			figure(B.fig[:number])
			clf()
		end
		return Void
	end

	fig["canvas"]["mpl_connect"]("pick_event", handlerPick)
	fig["canvas"]["mpl_connect"]("key_press_event", handlerKey)

	@schedule while true
		try B.fig[:show]() catch ex return end
		try plotBifurcation(B) end
		sleep(1)
	end

	return B
end

#update changed plot lines
function plotBifurcation(B::BifPlot)
	figure(B.fig[:number])
	subplot(121)

	cacheNew = Dict()

	#TODO parallel, changing length of projection
	lock(lockProject)
	for branch in B.project.branches
		if !haskey(B.cache, branch.solutions)
			x = map(last, branch.solutions)
			y = reduce(map(transpose∘projection, branch.solutions)) do A,a
				sA,sa = size(A,2), size(a,2)
				s = max(sA,sa)
				return [ A[:,(0:s-1)%sA+1]; a[1,(0:s-1)%sa+1] ]
			end

			lines = plot(x, y, picker=5, color="k", marker=".", markersize=3)
			cacheNew[branch.solutions] = lines
			map(l -> (B.idx[l]=branch), lines)
		else
			cacheNew[branch.solutions] = B.cache[branch.solutions]
			delete!(B.cache, branch.solutions)
		end
	end

	S = B.project.activeSolution
	if S ≠ Void
		if !haskey(B.cache, S)
			cacheNew[S] = plotSolution(B, S)
		else
			cacheNew[S] = B.cache[S]
			delete!(B.cache, S)
		end
	end
	unlock(lockProject)

	#cache now contains only obsolete plots
	for (k,v) in B.cache map(l->l["remove"](), v) end
	B.cache = cacheNew

	PyPlot.draw()

	return Void
end


function plotSolution(B::BifPlot, V::Vector{Float64})
	figure(B.fig[:number])
	subplot(121)

	x = V[end]
	y = projection(V)
	p1 = scatter(fill(x, size(y)), y, color="r")

	subplot(122, projection="3d")
	t = linspace(0, 2pi, length(V)÷3 * 8)
	v = reduce(hcat, toTrajInterp(V, 3)(t))'
	p2 = plot(v[:,1], v[:,2], v[:,3], color="k")

	return [p1;p2]
end
