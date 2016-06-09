# abstract type representing homotopies
# for T<:Homotopy H(::T), J(::T), V(::T), C(::T) MUST be overridden.
#	H - R^(N+1)->R^N function
#	J - R^(N+1)->R^(Nx(N+1)) function
#	V - view function. solution->appropriate view
#	C - control function. ()->controls
#	S - toScalar function. solution->[scalars] for bifurcation view

# C is executed through  View - Single Solution View  o.ä.
# V is executed through clicking a solution

#TODO new system callback

#TODO name aptly
#TODO howto forward difference, etc?
#	whole S window is application dependent
#	what if continuation method does not need J?

abstract Homotopy

H(x::Homotopy) = error("H not defined for $(typeof(x))")
J(x::Homotopy) = error("J not defined for $(typeof(x))")
V(x::Homotopy) = error("V not defined for $(typeof(x))")
C(x::Homotopy) = error("C not defined for $(typeof(x))")
S(x::Homotopy) = error("S not defined for $(typeof(x))")

# type encapsulating solutions and their sequence
# this is bifurcation view territory
type branch
	S::Vector{Vector{Float64}}
end
