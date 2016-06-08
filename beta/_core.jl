# this file is the entry point for the program.
# whole program is structured MVC

# paradigms
#	the program itself only encapsulates continuation stuff (no galerkin etc.)
#	functionality usable as library

# M
#	shall be equivalent to a saved file (branches, type of system, ...)
#	the only way to change the model is through continuation steps, as the model consists only of solutions and their order.
#	model is a singleton object, existence is a prerequisite for the rest of the program

# V
#	is static display of defect and time domain of the selected solution on one and bifurcation view on the other.
#	however, might opt to include customization feature
#	time domain view not mandatory... can use only continuation

# C
#	realized through controls of a type of system



#TODO interface mock-up till fr! (interactivity but static data)

# Qs
#	interface C->M?
#		C generates solutions, commit to M
#		consider current branch WIP?
#			ensure safety of other branches
#			duplicate branch? undo?
#			commit finished branch?
#	initial finding



# B window
#	set limits
#	set perturbation?
#	start continuation

# T window
#	only necessary to view/manipulate single solution
#	can create new branches from here?

include("GUI.jl")

startGUI()
