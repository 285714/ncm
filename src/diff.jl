
function forwardDiff(f, ɛ = 1e-8)
	X -> begin
		fX = f(X);
		h  = length(fX)
		w  = length(X)
		J  = zeros(h, w)
		foldl((J,i) -> [J   ( f(X + ɛ .* ((1:w).==i)) - fX ) ./ ɛ], zeros(h,0), 1:w)
	end
end


function interpolDiff(f, n = 5, h₀ = 0.01)
	hs = (0:h₀/n:h₀)[2:end]
	# hs = h₀ .* cos( π .* (2 .* (1:n) + 1) ./ (4n + 2) )

	X -> begin
		fX = f(X)
		h  = length(fX)
		w  = length(X)
		J  = zeros(h, w)

		for i = 1:w
			ϵ = (1:w) .== i
			ϕ(h) = ( f(X + h .* ϵ) - f(X - h .* ϵ) ) ./ 2h
			neville = zeros(h, n)

			for k = 0:n
				for j = 1:n-k
					if k == 0
						neville[:,j] = ϕ(hs[j])
					else
						neville[:,j] = (( -hs[j+k].^2 .* neville[:,j] + hs[j].^2 .* neville[:,j+1] )
							./ ( hs[j].^2 - hs[j+k].^2 ))
					end
				end
			end

			J[:,i] = neville[:,1]
		end

		J
	end
end





using PyPlot

n = 10
g(x)  = exp(cos( (1:n) .* sum(x) ))
Dg(x) = begin s = (1:n) .* sum(x) ; ( -exp(cos(s)) .* sin(s) .* (1:n) ) * ones(1,n) end
x = rand(10)
Y = zeros(0,2)

for i = 0:18
	h = 1/10^i

	fDg = forwardDiff(g, h)
	iDg = interpolDiff(g, 5, h)

	J  = Dg(x)
	fJ = fDg(x)
	iJ = iDg(x)

	fR = norm(J - fJ)
	iR = norm(J - iJ)

	Y = [Y; fR iR]
end

ion()
plot(log(Y[:,1]), label="forwardDiff")
plot(log(Y[:,2]), label="interpolDiff")
legend(loc="upper right",fancybox="true")
xlabel("1/10^i")


