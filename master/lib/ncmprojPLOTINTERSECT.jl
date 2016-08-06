function plotintersect(S::Session, c; its=20000, h=.05, rkonly=false)
	local data
	rk4((t,x)->S.core.f(t,[x;c]), .0, 10rand(3), h, predCount(its), init=()->(data=Array{Float64}(0,3)), callback=(t,y,f)->(data = [data;y']))

	f = figure()
	axes(projection="3d")
	hold(true)
	plot(data[end÷4:end,1], data[end÷4:end,2], data[end÷4:end,3], color="k", alpha=rkonly?1.0:.3)
	f[:show]()
	PyPlot.draw()
	rkonly && return
	for B in S.P
		col = get(S.viz.color, B, "k")
		for (a,b) in zip(B,B[2:end])
			if sign(a[end] - c) ≠ sign(b[end] - c)
				# k = abs((a[end]-c)/(a[end]-b[end]))
				# x = k.*a.data + (1-k).*b.data
				x = abs(a[end]-c) < abs(b[end]-c) ? a.data : b.data
				Htmp(V) = S.core.H([V; c])
				Jtmp(V) = S.core.J([V; c])[:, 1:end-1]
				x = newton(Htmp, Jtmp, x[1:end-1], predCount(25) ∧ predEps(1e-8))
				x = [x;c]
				X,ω,ℵ = unwrap(x)
				m = size(X,1)-1
				f = interp(X)
				T = linspace(0,2pi,8m)
				F = f(T)
				plot(F[:,1], F[:,2], F[:,3], color=col, alpha=.7)
			end
		end
	end
	PyPlot.draw()
end
