#simplified concept
#	two components: continuation & system
#	no interface requirements: components call stuff and it either exists or it doesn't
#		exceptions: continuationExec, systemExec must exist
#		types only for datatypes (if at all), generalization is user responsibility
#		global Proj variable is accessed by both components
#			never simultaneously, no sync, no communication (in general, otherwise: user responsibility)
#		no contraints on what components are allowed to do: user responsibility

#	this file specifies the whole interface...

include("Model.jl")

addprocs(1)

push!(LOAD_PATH, "$(pwd())/lib")
using Gtk.ShortNames

#(do menu stuff, other global GUI stuff, saving, ...)

#create empty project: vector of branches of solutions
global project = Project(Branch[])
#(or load project data)

#select continuation method, select system
# include(open_dialog("Select continuation method.", Gtk.GtkNullContainer(), ASCIIString[]))
# include(open_dialog("Select system.", Gtk.GtkNullContainer(), ASCIIString[]))

#DEBUG prevent loading screens
include("continuation/PC.jl")
include("system/roessler.jl")

# start
continuationExec()
systemExec()


#TODO-2DAY:
#	choose branch ✓
#	select solution ✓
#	add branch ✓
#	modify solution
#	delete branches
#	plot thread... not working... forget that... pyplot needs thread 1... Gtk might work...
#	singular branches...

#	abstraction: plotting, galerkin, continuation, project (!)
#	howto associate additional info with solutions/branches? transparent, general...
#		projects must be method agnostic...
#		plots, projection, stepsize, ...
#		encapsulate to simulate inheritance? serializability?
