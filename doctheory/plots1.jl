const a4width = 8.2677165354
const a4height = 11.6929133858

for b in ses.P
	for s in b
		if s[end] > 25; unshift!(b)
		else break end
	end

	for s in b[end:-1:1]
		if s[end] > 25; pop!(b)
		else break end
	end
end

f = figure(ses.viz.figBif[:number])
f[:set_figwidth](a4height)
f[:set_figheight](a4width)
f[:set_dpi](200)
PyPlot.draw()
savefig("roessleroverview.svg", format="svg")


map(println, keys(f))

# cut
include("$MASTER/lib/ncmprojPLOTINTERSECT.jl")
plotintersect(ses, 20, its=80000, h=.01)

f = gcf()
f[:set_figwidth](a4height)
f[:set_figheight](a4width)
f[:set_dpi](200)
PyPlot.draw()
savefig("roesslercut20.svg", format="svg")






# lorenz

f = figure(ses.viz.figBif[:number])
f[:set_figwidth](a4height)
f[:set_figheight](a4width)
f[:set_dpi](50)
PyPlot.draw()
savefig("lorenzoverview.svg", format="svg")


map(println, keys(f))

# cut
include("$MASTER/lib/ncmprojPLOTINTERSECT.jl")
plotintersect(ses, 80, its=80000, h=.005)

f = gcf()
f[:set_figwidth](a4height)
f[:set_figheight](a4width)
f[:set_dpi](200)
PyPlot.draw()
savefig("lorenzcut80.svg", format="svg")
