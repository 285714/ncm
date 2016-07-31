"""
todo
"""
type Solution
	data::Vector{Float64}
	parent # ::Branch # circular type dependency...
	Solution(d::Vector{Float64}) = new(d,nothing)
	Solution(d::Vector{Float64}, p) = new(d,p)
	Solution(s::Solution, p) = new(deepcopy(s.data), p)
end
@relay(Solution, :data)

Base.show(io::Base.IO, P::Solution) = write(io, "***SOLUTION***")
