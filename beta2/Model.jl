# type Solution
# 	raw::Vector{Float64}
# end

type Branch
	solutions::Vector{Vector{Float64}}	# ordered
	D::Dict #general storage... used for branch specific data
	Branch(solutions::Vector{Vector{Float64}}) = new(solutions, Dict())
	Branch(solution::Vector{Float64}) = new(Vector{Float64}[solution], Dict())
end

type Project
	branches::Vector{Branch}	# unordered
end
