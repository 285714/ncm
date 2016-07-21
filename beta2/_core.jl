#simplified concept
#	two components: continuation & system
#	no interface requirements: components call stuff and it either exists or it doesn't
#		exceptions: continuationExec, systemExec must exist
#		types only for datatypes (if at all), generalization is user responsibility
#		global Proj variable is accessed by both components
#			never simultaneously, no sync, no communication (in general, otherwise: user responsibility)
#		no contraints on what components are allowed to do: user responsibility

#	this file specifies the whole interface...

#TODO ifft for trigoninterp

addprocs(1)

push!(LOAD_PATH, "$(pwd())/lib")
using Gtk.ShortNames
@everywhere using mbNewton, mbPred, mbUtil, mbInterpolate, mbObserve
@everywhere include("Model.jl")
include("lib/ncmprojMKCONTROLGRID.jl")

#(do menu stuff, other global GUI stuff, saving, ...)

#create empty project: vector of branches of solutions
global project = Project()
global lockProject = ReentrantLock()

save(filename; overwrite=false) = if !isfile("save/$(filename)") || overwrite
	open("save/$(filename)", "w") do f serialize(f, project) end
else error("File already exists. Use  overwrite=true  .") end
load(filename) = open(deserialize, "save/$(filename)")

#select continuation method, select system
# include(open_dialog("Select continuation method.", Gtk.GtkNullContainer(), ASCIIString[]))
# include(open_dialog("Select system.", Gtk.GtkNullContainer(), ASCIIString[]))

@everywhere include("system/lorenz.jl")
# @everywhere include("system/roessler.jl")

include("lib/ncmprojBIFPLOT.jl") #needs thread 1
@schedule B = BifPlot(project)

@everywhere include("continuation/PC.jl")
@schedule pcGUI()

@everywhere include("lib/ncmprojGALERKIN.jl")
@schedule galerkinGUI()


#TODO-2DAY:
#	abstraction: plotting ✓, galerkin, continuation, project (! - for access, sync, sanity, ...)
#		lock
#		project: all operations
#	howto associate additional info with solutions/branches? transparent, general...
#		projects must be method agnostic...
#		plots, projection, stepsize, ...
#		encapsulate to simulate inheritance? serializability?

# (perturbation)
# bug: stepsize is associated with branch, not correct on deletion... reset, make changeable, associate with solution
# bug: apparently activeSoltion read from wrong (global) project...

# poincare ✓
# improv RK	?
# perturbation ✓
# lorenz

#TODO keys: delete from here, delete till here
#TODO h box, fix
#TODO broyden
#TODO mark current branch


# plot:
#	mark curr branch ✓
#	color ✓
#	lock branch (make unselectable, how to unlock?, how switch branches efficiently?)
#	rebuild completely as observable


# winkel pred ✓
# speichern -> encapsulate galerkin, make serializable
# lorenz
# roessler bif

#mark perturb branches
#auto new branch
#singleton perturb branch?
#new panel?
#special perturb mode?
#change perturb style... make factor for rand

#broyden pc, keep J per branch


#=
c = 350.0
m = 64
cyc, P, ω = findCyclePoincare((t,v)->f(t,[v;c]), rand(3),
	nIntersections=120, maxCycles=30, sampleSize=m,
	transientIterations=4000, transientStepSize=.01,
	steadyStateStepSize=.01)

C = resample(rfft(cyc, [1]), m)
C = [vec(vcat(real(C), imag(C[2:end,:]))); ω]

D = [C;c]
setActiveSolution(project, D)
d = [100*rand(size(D,1)-2); 0; 0]
@time newton(H, J, D, predCount(6), callback=(v,H,J)->println(norm(H)))
@time newton(Hroessler, Jroessler, D, predCount(6), callback=(v,H,J)->println(norm(H)))
@time newton(Hroessler, forwardDifference(Hroessler), D, predCount(6), callback=(v,H,J)->println(norm(H)))

s(x,r=Inf) = begin matshow(clamp(x,-r,r)); colorbar() end
=#
