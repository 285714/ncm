# relay execution of certain functions on type  T  to one of its fields  f  .
# kind of a workaround for missing inheritance
#TODO rethink this... might fuck up interface e.g. in conjuction with observer
macro relay(T′, f′)
	T,f = eval(T′), eval(f′)	# use esc instead?
	@assert f in fieldnames(T)
	U = fieldtype(T, f)

	return quote
		Base.setindex!(x::$T, v...) 	= setindex!(x.$f, v...)
		Base.getindex(x::$T, v...) 		= getindex(x.$f, v...)
		Base.start(x::$T) 				= start(x.$f)
		Base.next(x::$T, state) 		= next(x.$f, state)
		Base.done(x::$T, state) 		= done(x.$f, state)
		Base.convert(::Type{$U}, x::$T)	= x.$f
		Base.push!(x::$T, items...)		= push!(x.$f, items...)
		Base.pop!(x::$T, items...)		= pop!(x.$f, items...)
		Base.unshift!(x::$T, items...)	= unshift!(x.$f, items...)
		Base.shift!(x::$T, items...)	= shift!(x.$f, items...)
		Base.length(x::$T)				= length(x.$f)
		Base.endof(x::$T)				= endof(x.$f)
		return Void
	end
end

#operations
#	find branch of solution
#	add solution to branch
#	del solution from branch
# 	add/del branch
#	iterate over solutions of branch
#	iterate over branches
#	find other solutions of branch (find branch, iterate over branch)
#	manipulate solution? consider different solution? homo-/heterogeneous branches?
#	associate information with solutions, branches... -> responsibility of other subsystems
#		facilate through observer

## SOLUTION ##

type Solution
	data::Vector{Float64}
	parent # ::Branch # circular type dependency...
	Solution(d::Vector{Float64}, p) = new(d,p)
	Solution(s::Solution, p) = new(deepcopy(s.data, p))
end
@relay(Solution, :data)

## BRANCH ##

#TODO fix contructors

type Branch
	solutions::Vector{Solution}	# ordered
	parent # ::Project
	Branch(V::Vector{Solution}) = new(V, Void)
	Branch(P) = Branch(Solution[], P) #TODO runtime typeassert bco circ type?
	Branch(S,P) = new(S,P)
end
@relay(Branch, :solutions)

function Base.push!(B::Branch, sols...)
	P = B.parent
	for s in sols
		push!(B.solutions, s)
		mbObserve.notify(P, :pushSolution, B, s)
	end
	return P
end

function Base.unshift!(B::Branch, sols...)
	P = B.parent
	for s in sols
		push!(B.solutions, s)
		mbObserve.notify(P, :unshiftSolution, B, s)
	end
	return P
end

function Base.pop!(B::Branch)
	P = B.parent
	s = pop!(B.solutions)
	mbObserve.notify(P, :popSolution, B, s)
	return s
end

function Base.shift!(B::Branch)
	P = B.parent
	s = shift!(B.solutions)
	mbObserve.notify(P, :shiftSolution, B, s)
	return s
end

## PROJECT ##

type Project
	branches::Vector{Branch}	# unordered
	activeSolution::Union{Solution, Vector{Float64}, Type{Void}}
	O::Observer
	lock::ReentrantLock
	Project() = new(Branch[], Void, Observer(), ReentrantLock())
	Project(B::Vector{Branch}) = new(B, Void, Observer(), ReentrantLock())
end
@relay(Project, :branches)

#TODO lots of functions missing... add? pop e.g. really unlikely op... -> add on demand

Base.lock(P::Project) = lock(P.lock)
Base.unlock(P::Project) = unlock(P.lock)

mbObserve.observe(cb::Function, P::Project, s) = observe(cb, P.O, s)
mbObserve.notify(P::Project, s, v...) = mbObserve.notify(P.O, s, v...)

function Base.push!(P::Project, branches...)
	for b in branches
		push!(P.branches, b)
		mbObserve.notify(P, :pushBranch, b)
	end
	return P
end

function setActiveSolution(P::Project, V::Vector{Float64})
	P.activeSolution = V
	mbObserve.notify(P, :activeSolutionChanged)
	return Void
end

function setActiveSolution(P::Project, S::Solution)
	B = S.parent
	P = B.parent
	P.activeSolution = S
	mbObserve.notify(P, :activeSolutionChanged)
	return Void
end

function findSolution(P::Project, S::Vector{Float64})
	for i in 1:length(P.branches)
		j = findfirst(x -> x == S, P.branches[i].solutions)
		j > 0 && return i,j
	end
	return 0,0
end
