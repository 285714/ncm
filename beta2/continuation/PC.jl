using Gtk.ShortNames
using matsboNWTN
include("../lib/ncmprojMKCONTROLGRID.jl")
include("../lib/ncmprojLOGGER.jl")

function pcGUI()
	win = @Window("Continuation Controls", 256, 256, false)
	# setproperty!(win, :resizable, false)

	g = @Grid()
	push!(win, g)

	# global L = Logger() #TODO encapsulate
	# g[1:2,1] = L.w

	dataPC, gridPC = mkControlGrid([
		("ϵ", Float64, 1e-6, 1e-2, 1e-8),
		("κ", Float64, .0, 1.0, 1e-3),
		("δ", Float64, .0, 10.0, 1e-2),
		("α", Float64, .0, 90.0, 1e-1)
	])
	g[1:2,2] = gridPC

	buttonStep = @Button("Step")
	g[1,3] = buttonStep
	signal_connect(buttonStep, "clicked") do w
		@schedule begin
			project.activeSolution == Void && return Void
			setproperty!(w, :sensitive, false)
			goSingleStep(dataPC)
			setproperty!(w, :sensitive, true)
		end
	end

	buttonRun = @ToggleButton("Run")
	g[2,3] = buttonRun
	signal_connect(buttonRun, "toggled") do w
		if getproperty(w, :active, Bool)
			@schedule while getproperty(w, :active, Bool)
				goSingleStep(dataPC)
				yield()
			end
		end
	end

	showall(win)
end

global h = Dict() #TODO encapsulate, keep per solution and branch
function goSingleStep(D)
	lock(lockProject)
	i,j = @fetch findSolution(project, project.activeSolution) #TODO hint parent
	V = project.activeSolution

	local branch,htmp
	if i ≠ 0 && j == length(project.branches[i].solutions)
		branch = project.branches[i]
		htmp = get(h, (branch, true), 1)
	elseif i ≠ 0 && j == 1
		branch = project.branches[i]
		htmp = get(h, (branch, false), -1)
	else #solution not at start/end of branch, not in branch: create new branch, continue in pos dir
		htmp = 1
		V = deepcopy(V)
		branch = Branch(V)
		push!(project.branches, branch)
		println("new branch!") #TODO remove
	end

	W,htmp = @fetch continuationStep(H, J, V, htmp, D["ϵ"], D["κ"], D["δ"], D["α"])

	# write(L, htmp)

	# write back
	dir = htmp > 0
	(dir?push!:unshift!)(branch.solutions, W)
	project.activeSolution = W
	h[(branch, dir)] = htmp

	unlock(lockProject)
	return Void
end

# untested
@everywhere function continuationStep(H, J, V, h, ϵ, κ, δ, α)
	# tangent-vector for matrix A
	function tang(A)
		t = nullspace(A)[:,1] #[A; ones(1,size(A,2))] \ [zeros(size(A,1)); 1]
		flipsign(t / norm(t), det([A; t']))
	end

	# newton corrector step
	function Δ(v)
		Hv, Jv = H(v), J(v)
		Hv, Jv, Jv \ Hv
	end

	local v₁, Hv₁
	while true
		v₀ = V + h * tang(J(V))

		# step size control
		Hv₀, Jv₀, Δv₀ = Δ(v₀)
		v₁ = v₀ - Δv₀
		Hv₁, Jv₁, Δv₁ = Δ(v₁)
		v₁ = v₁ - Δv₁
		κ′ = norm(Δv₁) / norm(Δv₀)
		δ′ = norm(Δv₀)
		α′ = acos( dot(tang(Jv₁), tang(Jv₀)) )

		f = clamp(max(sqrt(κ′/κ), sqrt(δ′/δ), α′/α), .5, 2.0)
		h /= f

		.5 < f < 2.0 && break
	end

	v₁ = newton(H, J, v₁, predEps(ϵ))

	return v₁, h
end
