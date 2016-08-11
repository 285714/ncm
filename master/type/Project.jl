"""
Contains all found `Branch`es
"""
type Project
	branches::Vector{Branch}	# unordered
	activeSolution::Union{Solution, Vector{Float64}, Type{Void}}
	O::Observer
	# lock::ReentrantLock
	Project() = new(Branch[], Void, Observer())
	Project(B::Vector{Branch}) = new(B, Void, Observer())
end
@relay(Project, :branches)

#prevent enormous amounts of data being printed
Base.show(io::Base.IO, P::Project) = write(io, "***PROJECT***")

#TODO lots of functions missing... add? pop e.g. really unlikely op... -> add on demand

# Base.lock(P::Project) = lock(P.lock)
# Base.unlock(P::Project) = unlock(P.lock)

mbObserve.observe(cb::Function, P::Project, s) = observe(cb, P.O, s)
mbObserve.notify(P::Project, s, v...) = mbObserve.notify(P.O, s, v...)

function Base.push!(P::Project, branches...)
	for b in branches
		push!(P.branches, b)
		mbObserve.notify(P, :pushBranch, b)
	end
	return P
end

function Base.deleteat!(P::Project, i::Int)
	B = P[i]
	P.activeSolution in B && setActiveSolution(P, P.activeSolution.data)
	deleteat!(P.branches, i)
	mbObserve.notify(P, :delBranch, B)
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
