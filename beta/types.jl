# abstract type representing homotopies
# for T<:Homotopy H(::T), J(::T), V(::T), C(::T) MUST be overridden.
#	H - R^(N+1)->R^N function
#	J - R^(N+1)->R^(Nx(N+1)) function
#	V - view function. solution->appropriate view
#	C - control function. ()->controlss

#TODO name aptly
#TODO howto forward difference, etc?
#TODO combine C,V - both reside on the same pane and are implementation specific
#	whole S window is application dependent

abstract Homotopy

H(x::Homotopy) = error("H not defined for $(typeof(x))")
J(x::Homotopy) = error("J not defined for $(typeof(x))")
V(x::Homotopy) = error("V not defined for $(typeof(x))")
C(x::Homotopy) = error("C not defined for $(typeof(x))")
