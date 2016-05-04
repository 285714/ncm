include("lib/matsboPRED.jl")

p1 = matsboPRED.predEps(.1)
p2 = matsboPRED.predCount(5)

∧(p...) = function pred(v...; init::Bool=false)
	return foldl( (A,a) -> A&&a, map( f->f(v..., init=init), p ) )
end

p = p1 ∧ p2

p(5)
