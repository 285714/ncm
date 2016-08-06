#simplified concept
#	two components: continuation & system
#	no interface requirements: components call stuff and it either exists or it doesn't
#		exceptions: continuationExec, systemExec must exist
#		types only for datatypes (if at all), generalization is user responsibility
#		global Proj variable is accessed by both components
#			never simultaneously, no sync, no communication (in general, otherwise: user responsibility)
#		no contraints on what components are allowed to do: user responsibility

#	this file specifies the whole interface...

#TODO overarching  Session  type, encapsulating Project, SystemCore, Continuation, Visualization
#	make activeSolution a property of session? removes dependency between SystemCore and Project

#TODO ifft for trigoninterp

# addprocs(1)

const MASTER = splitdir(@__FILE__())[1]
push!(LOAD_PATH, "$MASTER/lib")
using Gtk.ShortNames
using mbNewton, mbPred, mbUtil, mbInterpolate, mbObserve
map(include, [
	"lib/ncmprojMKCONTROLGRID.jl"; "lib/ncmprojRelay.jl";
	"type/Solution.jl"; "type/Branch.jl"; "type/Project.jl";
	"type/SystemCore.jl"; "type/ContinuationMethod.jl"; "type/Visualization.jl";
	"type/Session.jl";
	"type/PC.jl"; "type/Galerkin.jl"; "type/GalerkinViz.jl"
	"lib/ncmprojSAVEHELPER.jl"
	])

#(do menu stuff, other global GUI stuff, saving, ...)

#select continuation method, select system
# include(open_dialog("Select continuation method.", Gtk.GtkNullContainer(), ASCIIString[]))
# include(open_dialog("Select system.", Gtk.GtkNullContainer(), ASCIIString[]))

# @everywhere include("system/lorenz.jl")
# sesL = create(lorenz, lorenz′, lorenzProjection)
# sesL = load("lorenz350nodoublings", lorenz, lorenz′, lorenzProjection)
# save("lorenz", sesL, overwrite=true)

# @everywhere include("system/roessler.jl")
# ses = create(roessler, Hroessler, Jroessler, roesslerProjection)
# ses = create(roessler, roessler′, roesslerProjection)
# ses = load("...", roessler, Hroessler, Jroessler, roesslerProjection)
# ses = load("roessler4612126", roessler, Hroessler, Jroessler, roesslerProjection)
# save("roessler4612126", ses, overwrite=true)

#=
using mbRK
global data
rk4((t,x)->ses.core.f(t,[x;100.5]), .0, rand(3), .01, predCount(10000), init=()->(global data=Array{Float64}(0,3)), callback=(t,y,f)->(global data = [data;y']))

ion()
f = figure()
axes(projection="3d")
plot(data[end÷2:end,1], data[end÷2:end,2], data[end÷2:end,3])
=#


#=
using ProfileView
Profile.clear()
@profile for i in 1:10
	step(ses.cont)
end
ProfileView.view()
=#

# homokline orbits
# thompsen steward
# index fixpoint -> perturbation
# proj T
