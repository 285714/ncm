#TODO keep max/min, prevent relim

using Munkres
using PyCall
pygui(:tk) #prevent conflict with gtk
using PyPlot
# remove all keybindings
D = PyCall.PyDict(matplotlib["rcParams"])
for k in keys(D); contains(k, "keymap") && (D[k] = Any[]) end
ioff()

#####

# this visualization module expects functions specific to Galerkin (projection)

"""
todo
"""
type GalerkinViz <: Visualization
	parent::Session
	figBif
	figSol
	idxLineBranch::Dict
	idxBranchLines::Dict
	color::Dict
	activeSolutionMark
end

function Base.show(V::GalerkinViz)
	figure(V.figBif[:number])
	clf()
	PyPlot.draw()
	V.idxBranchLines = Dict{Any,Any}()
	V.idxLineBranch = Dict{Any,Any}()
	map(b -> try pushBranch(V,b) end, V.parent.P)
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

	V = GalerkinViz(parent, figBif, figSol, Dict{Any,Any}(), Dict{Any,Any}(), Dict{Any,Any}(), Void)

	figBif["canvas"]["mpl_connect"]("pick_event", ev -> handlerPick(V,ev))
	figBif["canvas"]["mpl_connect"]("key_press_event", ev -> @schedule handlerKey(V,ev))

	#connect to proj events
	P = parent.P

	observe(P, :pushBranch) do b
		pushBranch(V, b) end

	observe(P, :activeSolutionChanged) do
		s = V.parent.P.activeSolution #TODO no direct access...
		isa(s, Solution) && plotSolution(V, s.data)
		isa(s, Vector{Float64}) && plotSolution(V, s)
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



@noinline function handlerPick(V::GalerkinViz, ev) b = get(V.idxLineBranch, hash(ev[:artist]), nothing)
	b == nothing && return

	try
		i = ev[:ind][1]+1
		setActiveSolution(V.parent.P, b[i])
	catch e
		println(e)
	end
	return Void
end

@noinline function handlerKey(V::GalerkinViz, ev)
	# println(ev[:key]) #debug

	figure(V.figBif[:number])

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
		# figure(V.figBif[:number])
		# clf()
		# PyPlot.draw()
		# sleep(.1)
		# V.idxBranchLines = Dict{Any,Any}()
		# V.idxLineBranch = Dict{Any,Any}()
		# map(b -> pushBranch(V,b), V.parent.P)
		show(V)
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
		for L in values(V.idxBranchLines)
			for l in L
				l[:set_picker](0)
				l[:set_marker]("None")
			end
		end
	elseif ev[:key] == "ctrl+L"
		for L in values(V.idxBranchLines)
			for l in L
				l[:set_picker](5)
				l[:set_marker](".")
			end
		end
	elseif ev[:key] == "ctrl+h"
		!isa(V.parent.P.activeSolution, Solution) && return Void
		b = V.parent.P.activeSolution.parent
		for l in V.idxBranchLines[b]
			l[:set_picker](0)
			l[:set_visible](false)
		end
	elseif ev[:key] == "ctrl+H"
		for L in values(V.idxBranchLines)
			for l in L
				l[:set_picker](5)
				l[:set_visible](true)
			end
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
				V.color[V.parent.P.activeSolution.parent] = colMap[ev[:key]]
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
	elseif ev[:key] == "d"
		s = V.parent.P.activeSolution
		isa(s, Solution) && plotSolution(V, s.data; split=false)
		isa(s, Vector{Float64}) && plotSolution(V, s; split=false)
	elseif ev[:key] == "D"
		s = V.parent.P.activeSolution
		isa(s, Solution) && plotSolution(V, s.data; split=true)
		isa(s, Vector{Float64}) && plotSolution(V, s; split=true)
	end

	PyPlot.draw()
	return Void
end


function pushBranch(V::GalerkinViz, b::Branch)
	isempty(b) && error("empty branch")

	x = map(last, b)

	function reducer(A,a)
		sA,sa = size(A,2), size(a,2)

		if sA ≠ sa
			warn("in-branch change in number of projected points")
			tA,ta = vec(A[end,:]), vec(a)

			# find repetitions of elements of A s.t. norm(A[is]-B) is minimal
			function align(A,B)
				m,n = length(A), length(B)
				@assert m < n
				d(a,b) = norm(a-b)
				M = fill(Inf,m,n)
				M[1,1] = d(A[1],B[1])
				for j in 2:n-m+1 M[1,j] = M[1,j-1] + d(A[1],B[j]) end
				for i in 2:m, j in i:n-m+i
					M[i,j] = min(M[i-1,j-1], M[i,j-1]) + d(A[i],B[j])
				end
				is = [m]
				for j in n-1:-1:1
					i = is[1]
					unshift!(is, i ≠ 1 ? M[i-1,j] ≤ M[i,j] ? i-1 : i : 1)
				end
				return is
			end

			if sA>sa
				a = ta[align(ta,tA)]'
			else
				A = A[:,align(tA,ta)]
			end
		end

		return vcat(A,a)
	end

	function mapper(x)
		tmp = projection(V.parent.core, x)
		length(tmp) == 0 && warn("empty projection")
		return tmp'
	end

	y = mapreduce(mapper, reducer, b)

	figure(V.figBif[:number])
	hold(true)
	lines = plot(x, y, picker=5, color=get(V.color, b, "k"), marker=".", markersize=3)
	PyPlot.draw()

	map(l -> (V.idxLineBranch[hash(l)]=b), lines)
	V.idxBranchLines[b] = lines

	return lines
end


function addToBranch(op, V::GalerkinViz, b::Branch, s::Solution)
	!haskey(V.idxBranchLines, b) && error("branch missing")
	lines = V.idxBranchLines[b]

	x = last(s) #TODO extend relay?
	Y = projection(V.parent.core, s)

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
function plotSolution(V::GalerkinViz, v::Vector{Float64}; split=false)
	figure(V.figBif[:number])
	hold(true)
	try V.activeSolutionMark[:remove]() end
	x = v[end]
	y = projection(V.parent.core, v)
	V.activeSolutionMark = scatter(fill(x, size(y)), y; facecolors="none", edgecolors="k", marker="o")
	PyPlot.draw()


	figure(V.figSol[:number])
	cla()
	gca(projection="3d")
	axis("equal")
	hold(false)
	m = length(v)÷6
	t = linspace(0, 2pi, 2m*4)
	w = interp(unwrap(v)[1])(t)
	split && (w[:,3] += linspace(0, maximum(w[:,3])-minimum(w[:,3]), length(w[:,3])))
	plot(w[:,1], w[:,2], w[:,3], color="k")

	hold(true)
	try cbPlotSolution(v)
	catch e; println(e) end

	PyPlot.draw()

	return Void
end


function delBranch(V::GalerkinViz, b::Branch)
	!haskey(V.idxBranchLines, b) && return
	lines = V.idxBranchLines[b]
	delete!(V.idxBranchLines, b)
	map(lines) do l
		l[:remove]()
		delete!(V.idxLineBranch, hash(l))
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
		t = linspace(0, 2pi, length(v)÷3 * 6)
		w = interp(unwrap(v)[1])(t)
		plot(w[:,1], w[:,2], w[:,3], color="k", alpha=.1)
	end

	v = B[1]
	t = linspace(0, 2pi, length(v)÷3 * 6)
	w = interp(unwrap(v)[1])(t)
	plot(w[:,1], w[:,2], w[:,3], color="y")

	v = B[end]
	t = linspace(0, 2pi, length(v)÷3 * 6)
	w = interp(unwrap(v)[1])(t)
	plot(w[:,1], w[:,2], w[:,3], color="r")

	PyPlot.draw()
	return
end











#=
@noinline function hungarian(M)
	#make square
	addrows, addcols = 0,0
	if size(M,2) < size(M,1)
		addcols = size(M,1)-size(M,2)
		M = hcat(M, fill(maximum(M), size(M,1), addcols))
	elseif size(M,1) < size(M,2)
		addrows = size(M,2)-size(M,1)
		M = vcat(M, fill(maximum(M), addrows, size(M,2)))
	end

	L = size(M,1)

	#subtract row/col minima
	mapslices(c->c-minimum(c), M, [1])
	mapslices(r->r-minimum(r), M, [2])

	while true
		println(0)

		#count zeros
		is,js = zeros(L),zeros(L)
		zs = 0
		for i in 1:L, j in 1:L
			if M[i,j] == 0
				is[i] += 1
				js[j] += 1
				zs += 1
			end
		end

		#sorted idx
		seqi, seqj = sortperm(is, rev=true), sortperm(js, rev=true)

		println(1)

		#mapping with definite candidates set
		# mapping = zeros(Int64, L)
		# for k in L:-1:1
		# 	is[seqi[k]] == 1 && (mapping[seqi[k]] = findfirst(M[seqi[k],:], 0))
		# 	js[seqj[k]] == 1 && (mapping[findfirst(M[seqj[k],:]], 0)] = k)
		# 	is[seqi[k]] > 1 && js[seqj[k]] > 1 && break
		# end
		#
		# #check for solution
		# idx = find(x->x==0, mapping)
		# opts = [ find(x->x==0, M[i,:]) for i in idx]
		# K = length(opts)
		# count = zeros(Int64, K)
		# while true
		# 	#construct candidate solution, check
		# 	candidate = map((o,c)->o[c+1], opts, count)
		# 	flagBadCandidate = false
		# 	for c in 1:K
		# 		for d in c+1:K
		# 			candidate[c] == candidate[d] && (flagBadCandidate = true) && break
		# 		end
		# 		flagBadCandidate && break
		# 	end
		#
		# 	#good solution found
		# 	if !flagBadCandidate
		# 		mapping[idx] = candidate
		# 		return mapping
		# 	end
		#
		# 	#prepare next candidate
		# 	k = length(opts)
		# 	while true
		# 		count[k] = (count[k] + 1) % length(opts[k])
		# 		count[k] ≠ 0 && break
		# 		k -= 1
		# 		k == 0 && break
		# 	end
		# 	k == 0 && break
		# end

		#construct lines through zeros
		lis, ljs = [], []
		while zs > 0
			i,j = L,L
			while is[seqi[i]] == 1
				push!(lis, seqi[i])
				push!(lis, seqi[i])
			end
			if is[i] ≥ js[j]
				zs -= is[i]
				is[i] = 0
				seqi = circshift(seqi, [-1])
				push!(lis, i)
				for j in L:-1:1
					M[i,seqj[j]] ≠ 0 && continue
					js[seqj[j]] -= 1
					while j < L && js[seqj[j]] < js[seqj[j+1]]
						seqj[j],seqj[j+1] = seqj[j+1],seqj[j]
						j += 1
					end
				end
			else
				zs -= js[j]
				js[j] = 0
				seqj = circshift(seqj, [-1])
				push!(ljs, j)
				for i in L:-1:1
					M[seqi[i],j] ≠ 0 && continue
					is[seqi[i]] -= 1
					while i < L && is[seqi[i]] < is[seqi[i+1]]
						seqi[i],seqi[i+1] = seqi[i+1],seqi[i]
						i += 1
					end
				end
			end
		end

		println(2, " length(lis)+length(ljs) > L : ", length(lis)+length(ljs) > L)

		if length(lis)+length(ljs) == L
			println(ceil(Integer, M))
			opts = [ find(x->x==0, M[i,:]) for i in 1:L ]
			cnt = ones(Int64, L)
			for k in 1:prod(map(length, opts))
				mapping = map(getindex, opts, cnt)
				good = true
				for i in 1:L, j in i+1:L; mapping[i] == mapping[j] && (good = false; break) end
				good && return mapping[1:L-addcols]
				for i in 1:L
					cnt[i] = cnt[i] % length(opts[i]) + 1
					cnt[i] ≠ 1 && break
				end
			end
		end

		println(3)

		#prepare next round
		minunc = minimum(M[setdiff(1:L, lis), setdiff(1:L, ljs)])
		M[lis,:] += minunc
		M[:,ljs] += minunc
		M -= minunc

		println(4)
	end
end
hungarian(A,B) = hungarian([ norm(a-b) for a in A, b in B ])

hungarian(rand(5,5))
hungarian(rand(4,5))
hungarian(rand(5,4))

hungarian([1,4,4], [1,4])
=#
