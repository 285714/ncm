# type Solution
# 	data::Vector{Float64}
# end

type Branch
	solutions::Vector{Vector{Float64}}	# ordered
	h::Dict{Bool,Float64}
	Branch(solutions::Vector{Vector{Float64}}) = new(solutions, Dict(true=>1.0,false=>-1.0))
	Branch(solution::Vector{Float64}) = new(Vector{Float64}[solution], Dict(true=>1.0,false=>-1.0))
end

type Project
	branches::Vector{Branch}	# unordered
end
