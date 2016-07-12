using Gtk.ShortNames

function pcGUI()
	win = @Window("Continuation Controls", 256, 256, false)
	# setproperty!(win, :resizable, false)

	g = @Grid()
	push!(win, g)

	# global L = Logger() #TODO encapsulate
	# g[1:2,1] = L.w

	dataPC = Dict{AbstractString, Any}()
	gridPC = mkControlGrid(dataPC, [
		("ϵ", Float64, 1e-4, 1e-6:1e-8:1e-2),
		("κ", Float64, .4, .0:1e-3:1.0),
		("δ", Float64, .5, .0:1e-2:5.0),
		("α", Float64, 10.0, .0:.1:90.0),
		("Perturb", Float64, .0, -1e-3:1e-6:1e-3)
	])
	g[1:2,2] = gridPC

	buttonStep = @Button("Step")
	g[1,3] = buttonStep
	signal_connect(buttonStep, "clicked") do w
		@schedule begin
			project.activeSolution == Void && return Void
			setproperty!(w, :sensitive, false)
			lock(lockProject)
			try
				goSingleStep(dataPC)
			finally
				unlock(lockProject)
				setproperty!(w, :sensitive, true)
			end
		end
	end

	buttonRun = @ToggleButton("Run")
	g[2,3] = buttonRun
	signal_connect(buttonRun, "toggled") do w
		if getproperty(w, :active, Bool)
			@schedule while getproperty(w, :active, Bool)
				lock(lockProject)
				try goSingleStep(dataPC)
				finally unlock(lockProject) end
				sleep(.05)
			end
		end
	end

	showall(win)
end

global h = Dict() #TODO encapsulate, keep per solution and branch
function goSingleStep(D)
	V = project.activeSolution
	V == Void && return Void #TODO elaborate..

	local B
	if typeof(V) == Solution && (V == V.parent[1] || V == V.parent[end])
		B = V.parent
	else # internal, single or orphaned solution
		B = Branch(project)
		V = Solution(V, B)
		B.solutions = Solution[V]
		push!(project, B)
	end

	htmp = get!(h, V, 1.0)
	V==B[1] && length(B)>1 && (htmp=-htmp)

	Htmp = D["Perturb"]==0 ? H : V -> H(V) + D["Perturb"] #maybe prevent call overhead
	W,htmp = continuationStep(Htmp, J, V.data, htmp, D["ϵ"], D["κ"], D["δ"], D["α"]) #TODO fetch...

	W = Solution(W, B)
	h[W] = abs(htmp)

	(htmp>0?push!:unshift!)(B, W)
	setActiveSolution(project, W)

	return Void
end

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
		abs(h) ≤ 1e-8 && error("h≈0 !") #TODO remove
	end

	v₁ = newton(H, J, v₁, predEps(ϵ) ∧ predCount(25)) #TODO integrate max count
	return v₁, h
end
