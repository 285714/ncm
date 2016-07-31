#TODO deleting of solutions does not delete from h.

type PC <: ContinuationMethod
	parent::Session
	dataGUI::Dict{AbstractString, Any}
	h::Dict{Any,Float64}
	PC(ses::Session) = new(ses, Dict{AbstractString, Any}(), Dict{Any,Float64}())
end

function Base.show(C::PC)
	win = @Window("Continuation Controls", 256, 100, false)
	# setproperty!(win, :resizable, false)

	g = @Grid() #TODO obsolete
	push!(win, g)

	buttonStep = @Button("Step")
	buttonRun = @ToggleButton("Run")

	gridPC = mkControlGrid(C.dataGUI, Array[
		#[ ("ϵ", Float64, 1e-4, 1e-6:1e-8:1e-2), Void ],
		[ ("κ", Float64, .4, 1e-3:1e-3:1.0) ],
		[ ("δ", Float64, .5, 1e-2:1e-2:20.0) ],
		[ ("α", Float64, 10.0, .1:.1:40.0) ],
		[ ("Perturb", Bool, false); ("Strength", Float64, .0, -.1:1e-3:.1) ],
		[ buttonStep, buttonRun ]
	])
	g[1:2,2] = gridPC

	signal_connect(buttonStep, "clicked") do w
		@schedule begin
			C.parent.P.activeSolution == Void && return Void
			setproperty!(w, :sensitive, false)
			try
				step(C)
			finally
				setproperty!(w, :sensitive, true)
			end
		end
	end

	signal_connect(buttonRun, "toggled") do w
		if getproperty(w, :active, Bool)
			@schedule while getproperty(w, :active, Bool)
				try step(C)
				catch e
					setproperty!(w, :active, false)
					return
				end
				sleep(.2)
			end
		end
	end

	showall(win)
end

# can throw error, eg when h=0
function Base.step(C::PC)
	V = C.parent.P.activeSolution
	V == Void && return Void #TODO elaborate..

	local B, flagNewBranch = false
	if typeof(V) == Solution && (V == V.parent[1] || V == V.parent[end])
		B = V.parent
	else # internal, single or orphaned solution
		flagNewBranch = true
		B = Branch(C.parent.P)
		V = Solution(V, B)
		B.solutions = Solution[V]
	end

	htmp = get!(C.h, V, 1.0)
	V==B[1] && length(B)>1 && (htmp=-htmp)

	Htmp = C.dataGUI["Perturb"] ? V -> H(C.parent.core)(V) + C.dataGUI["Strength"] : H(C.parent.core)
	W,htmp = continuationStep(Htmp, J(C.parent.core), V.data, htmp, 1e-8, C.dataGUI["κ"], C.dataGUI["δ"], C.dataGUI["α"]) #fetch?

	W = Solution(W, B)
	C.h[W] = abs(htmp)

	flagNewBranch && push!(C.parent.P, B) #NOTE no branch pushed in case of error in continuationStep
	(htmp>0?push!:unshift!)(B, W)
	setActiveSolution(C.parent.P, W)

	return Void
end

# numerical core
function continuationStep(H, J, V, h, ϵ, κ, δ, α)
	# tangent-vector for matrix A
	function tang(A)
		t = [A; ones(1,size(A,2))] \ [zeros(size(A,1)); 1] # nullspace(A)[:,1] #takes up >90% processing time
		flipsign(t / norm(t), det([A; t']))
	end

	# newton corrector step
	function Δ(v)
		Hv, Jv = H(v), J(v)
		Hv, Jv, Jv \ Hv
	end

	local v₁, dir = nothing
	while true
		tJV = tang(J(V))
		v₀ = V + h * tJV

		# step size control
		Hv₀, Jv₀, Δv₀ = Δ(v₀)
		v₁ = v₀ - Δv₀

		Hv₁, Jv₁, Δv₁ = Δ(v₁)
		v₁ = v₁ - Δv₁
		κ′ = norm(Δv₁) / norm(Δv₀)
		δ′ = norm(Δv₀)
		α′ = acos( clamp(dot(tJV, tang(Jv₀)), -1, 1) ) #prevent domain exception through numerical error

		f = max(sqrt(κ′/κ), sqrt(δ′/δ), (α′/pi*180)/α)

		# prevent infinite loop: only allow changes of h in one direction
		dir == nothing && (dir = (f < 1.0))
		dir ≠ (f < 1.0) && break

		h /= clamp(f, .5, 2)
		.5 ≤ f ≤ 2.0 && break
		abs(h) ≤ 1e-8 && error("h≈0 !") #TODO remove?
	end

	its = 0
	v₁ = newton(H, J, v₁, predEps(ϵ) ∧ predCount(50); callback=(H,J,v)-> (its+=1)) #TODO integrate max count
	its ≥ 50 && warn("no convergence; error: ", norm(H(v₁)))
	return v₁, h
end
