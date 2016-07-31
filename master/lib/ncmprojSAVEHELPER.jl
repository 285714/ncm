simplifyProject(P::Project) = map(P) do branch
	map(branch) do solution
		solution.data
	end
end

#TODO actually convert(T, Project)
function complicateProject(V) #::Vector{Vector{Vector{Float64}}}
	P = Project()
	branches = map(V) do bData
		branch = Branch(P)
		solutions::Vector{Solution} = map(bData) do sData
			Solution(sData, branch)
		end
		branch.solutions = solutions
		return branch
	end
	P.branches = branches
	return P
end

"""
    create(homotopy, jacobian, projection)
"""
create(v...) = begin
	ses = Session()
	ses.P = Project()
	ses.core = Galerkin(ses, v...)
	ses.cont = PC(ses)
	ses.viz = GalerkinViz(ses)
	show(ses.cont); show(ses.core); show(ses.viz)
	return ses
end

"""
    save(filename, session[, overwrite])
"""
save(filename, S::Session; overwrite=false) = begin
	if !isfile("save/$(filename)") || overwrite; open("save/$(filename)", "w") do f serialize(f, simplifyProject(S.P)) end
	else error("File already exists. Use  overwrite=true  .")end
	return Void
end

#TODO restore non-serializable stuff (figures, observer)
"""
    load(filename, homotopy, jacobian, projection)
"""
load(filename, v...) = begin
	V = open(deserialize, "save/$(filename)")
	ses = Session()
	ses.P = complicateProject(V)
	ses.core = Galerkin(ses, v...)
	ses.cont = PC(ses)
	ses.viz = GalerkinViz(ses)
	show(ses.cont); show(ses.core); show(ses.viz)
	return ses
end
