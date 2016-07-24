type Branch
	solutions::Vector{Solution}	# ordered
	parent # ::Project
	Branch(V::Vector{Solution}) = new(V, Void)
	Branch(P) = Branch(Solution[], P) #TODO runtime typeassert bco circ type?
	Branch(S,P) = new(S,P)
end
@relay(Branch, :solutions)

Base.show(io::Base.IO, P::Branch) = write(io, "***BRANCH***")

#TODO IMPORTANT: handle case of activeSolution::Solution

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
		unshift!(B.solutions, s)
		mbObserve.notify(P, :unshiftSolution, B, s)
	end
	return P
end

function Base.pop!(B::Branch)
	P = B.parent
	s = pop!(B.solutions)
	s == P.activeSolution && (P.activeSolution=s.data)
	if length(B) == 1
		deleteat!(P, findfirst(P, B))
		B[1] == P.activeSolution && setActiveSolution(P, B[1].data)
		mbObserve.notify(P, :delBranch, B)
	else
		mbObserve.notify(P, :popSolution, B, s)
	end
	return s
end

function Base.shift!(B::Branch)
	P = B.parent
	s = shift!(B.solutions)
	s == P.activeSolution && (P.activeSolution=s.data)
	if length(B.solutions) == 1
		deleteat!(P, findfirst(P, B))
		B[1] == P.activeSolution && setActiveSolution(P, B[1].data)
		mbObserve.notify(P, :delBranch, B)
	else
		mbObserve.notify(P, :shiftSolution, B, s)
	end
	return s
end
