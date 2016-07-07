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

#(do menu stuff, other global GUI stuff, saving, ...)

#create empty project: vector of branches of solutions
global project = Project()
global lockProject = ReentrantLock()

save(filename) = open("save/$(filename)", "w") do f serialize(f, project) end
load(filename) = (global project = open(deserialize, "save/$(filename)"))

#select continuation method, select system
# include(open_dialog("Select continuation method.", Gtk.GtkNullContainer(), ASCIIString[]))
# include(open_dialog("Select system.", Gtk.GtkNullContainer(), ASCIIString[]))

#DEBUG prevent loading screens
@everywhere include("system/roessler.jl")

@everywhere include("continuation/PC.jl")
pcGUI()

@everywhere include("lib/ncmprojGALERKIN.jl")
galerkinGUI()

include("lib/ncmprojBIFPLOT.jl") #needs thread 1
B = BifPlot(project)


#TODO-2DAY:
#	choose branch ✓
#	select solution ✓
#	add branch ✓
#	modify solution
#	delete branches
#	plot thread... not working... forget that... pyplot needs thread 1... Gtk might work...
#	singular branches...

#	abstraction: plotting ✓, galerkin, continuation, project (! - for access, sync, sanity, ...)
#		lock
#		branch: dequeue only
#		project: all operations
#	howto associate additional info with solutions/branches? transparent, general...
#		projects must be method agnostic...
#		plots, projection, stepsize, ...
#		encapsulate to simulate inheritance? serializability?

# (perturbation)
# bug: stepsize is associated with branch, not correct on deletion... reset, make changeable, associate with solution
# bug: apparently activeSoltion read from wrong (global) project...

;

# poincare ✓
# improv RK	?
# perturbation ✓
# lorenz

#TODO keys: delete from here, delete till here
#TODO h box, fix
#TODO broyden
#TODO controlgrid custom dict
#TODO mark current branch


# plot:
#	mark curr branch ✓
#	color ✓
#	lock branch (make unselectable, how to unlock?, how switch branches efficiently?)
#	rebuild completely as observable
