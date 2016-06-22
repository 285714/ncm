# type Solution
# 	raw::Vector{Float64}
# end

type Branch
	solutions::Vector{Vector{Float64}}	# ordered
	Branch(solutions::Vector{Vector{Float64}}) = new(solutions)
	Branch(solution::Vector{Float64}) = new(Vector{Float64}[solution])
end

type Project
	branches::Vector{Branch}	# unordered
	activeSolution
	Project(B) = new(B, Void)
end

function findBranch(P::Project, S::Vector{Float64})
	for i in 1:length(P.branches)
		S in P.branches[i].solutions && return i
	end
	return 0
end
