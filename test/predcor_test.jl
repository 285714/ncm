include("../src/ncm.jl")
using PyPlot


# could be more complex (ContinuationMethod: v.v → v) and care about scope… aaand be in ncm.js
macro vecin(f)
	return :(v -> $f(v[1], v[2:end]))
end




n     = 10
g(x)  = exp(cos( (1:n) .* sum(x) ))
Dg(x) = begin s = (1:n) .* sum(x) ; ( -exp(cos(s)) .* sin(s) .* (1:n) ) * ones(1,n) end
F     = @vecin (λ,z) -> z - λ * g(z)
DF    = @vecin (λ,z) -> [-g(z)   eye(n) - λ * Dg(z)]
v₀    = zeros(n+1)

pc = PredictorCorrector()

#=let
	X = 0
	Y = 0
	cm(pc, F, DF, function(m) print("↣") ; X = [X; m.v[1]] ; Y = [Y; norm(m.v[2:end])] ; m.n < 1000 end , v₀)

	plot(X, Y)
end=#


t = Task(() -> cm(pc, F, DF, produce, v₀))

function runN(n :: Int)
	tic()
	i = 0
	for m in t
		println("$(m.n)> λ=$(m.v[1]),\t|v|=$(norm(m.v[2:end]))\t[κ=$(m.κ), δ=$(m.δ), α=$(m.α)]")
		(@vecin (λ,z) -> scatter(λ, norm(z)))(m.v)

		if (i+=1) >= n
			break
		end
	end
	print("$n steps in $(toc()) seconds")
end



