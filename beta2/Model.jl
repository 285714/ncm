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

function findSolution(P::Project, S::Vector{Float64})
	for i in 1:length(P.branches)
		j = findfirst(x -> x == S, P.branches[i].solutions)
		j > 0 && return i,j
	end
	return 0,0
end
