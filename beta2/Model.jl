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

findBranch(P::Project, S::Vector{Float64}) = for branch in P.branches; S in branch.solutions && return branch end
