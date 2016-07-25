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

push!(LOAD_PATH, "$(pwd())/lib")
using Gtk.ShortNames
@everywhere using mbNewton, mbPred, mbUtil, mbInterpolate, mbObserve
@everywhere map(include, [
	"lib/ncmprojMKCONTROLGRID.jl"; "lib/ncmprojRelay.jl";
	"type/Solution.jl"; "type/Branch.jl"; "type/Project.jl";
	"type/SystemCore.jl"; "type/ContinuationMethod.jl"; "type/Visualization.jl";
	"type/Session.jl";
	"type/PC.jl"; "type/Galerkin.jl"; "type/GalerkinViz.jl"
	])

#(do menu stuff, other global GUI stuff, saving, ...)

create() = begin
	ses = Session()
	ses.P = Project()
	ses.core = Galerkin(ses,lorenz,lorenz′)
	ses.cont = PC(ses)
	ses.viz = GalerkinViz(ses)
	show(ses.cont); show(ses.core); show(ses.viz)
	return ses
end

save(filename, S::Session; overwrite=false) = begin
	if !isfile("save/$(filename)") || overwrite; open("save/$(filename)", "w") do f serialize(f, S.P.branches) end
	else error("File already exists. Use  overwrite=true  .")end
	return Void
end

#TODO restore non-serializable stuff (figures, observer)
load(filename) = begin
	branches = open(deserialize, "save/$(filename)")
	ses = Session()
	ses.P = Project()
	for b in branches
		b.parent = ses.P
		length(b) < 2 && deleteat!(branches, findfirst(branches, b))
	end
	ses.P.branches = branches
	ses.core = Galerkin(ses,lorenz,lorenz′)
	ses.cont = PC(ses)
	ses.viz = GalerkinViz(ses)
	show(ses.cont); show(ses.core); show(ses.viz)
	return ses
end

#select continuation method, select system
# include(open_dialog("Select continuation method.", Gtk.GtkNullContainer(), ASCIIString[]))
# include(open_dialog("Select system.", Gtk.GtkNullContainer(), ASCIIString[]))

@everywhere include("system/lorenz.jl")
@everywhere include("system/roessler.jl")



#TODO keys: delete from here, delete till here
#TODO h box, fix
#TODO broyden
#TODO mark current branch

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



#operations
#	find branch of solution
#	add solution to branch
#	del solution from branch
# 	add/del branch
#	iterate over solutions of branch
#	iterate over branches
#	find other solutions of branch (find branch, iterate over branch)
#	manipulate solution? consider different solution? homo-/heterogeneous branches?
#	associate information with solutions, branches... -> responsibility of other subsystems
#		facilate through observer

# minimum constraint for h?

# for b in ses.P; length(b) < 2 && deleteat!(ses.P, findfirst(ses.P, b)) end
