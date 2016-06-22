# encapsulate plotting of bifurcation data
# provides a view on a project
# should be auto-updating for minimum interface

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


function BifPlot(project)
	fig = figure()
	fig["canvas"]["set_window_title"]("Bifurcation Plot")
	fig["show"]()

	#index for associating plotlines with solutions
	idx = Dict{Any,Any}()
	cache = Dict{Any,Any}()

	B = BifPlot(fig, project, idx, cache, Void)

	local scatterBranch
	function handlerBifurcationPick(ev)
		branch = B.idx[ev["artist"]]
		i = ev["ind"][1]+1

		plotSolution(B, branch.solutions[i])
	end

	fig["canvas"]["mpl_connect"]("pick_event", handlerBifurcationPick)

	@schedule while true #TODO prevent override
		plotBifurcation(B)
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
	for branch in B.project.branches
		if !haskey(B.cache, branch.solutions)
			x = map(last, branch.solutions)
			y = reduce(hcat, map(projection, branch.solutions))'
			lines = plot(x, y, picker=5, color="k", marker=".")

			# u = repmat(x[[1;end]], 1, size(y,2))
			# v = y[[1;end],:]
			# handles = scatter(u, v, color="k")

			cacheNew[branch] = lines #[lines; handles]
			map(l -> (B.idx[l]=branch), lines)
		else
			cacheNew[branch] = cache[branch]
			delete!(cache, branch)
		end
	end

	#cache now contains only obsolete plotlines
	for (k,v) in B.cache map(l->l["remove"](), v) end
	B.cache = cacheNew

	PyPlot.draw()

	return Void
end


function plotSolution(B::BifPlot, V)
	figure(B.fig[:number])
	subplot(121)

	try B.solutionMark["remove"]() end
	x = V[end]
	y = projection(V)
	B.solutionMark = scatter(fill(x, size(y)), y, color="r")

	subplot(122, projection="3d")
	t = linspace(0, 2pi, length(V)÷3 * 8)
	v = reduce(hcat, toTrajInterp(V, 3)(t))'
	plot(v[:,1], v[:,2], v[:,3], color="k")

	return Void
end
