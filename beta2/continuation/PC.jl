include("../lib/ncmprojMKCONTROLGRID.jl")

function continuationExec()
	continuationGUI()
end

function continuationGUI()
	win = @Window("Continuation Controls", 500, 100, false)

	g = @Grid()
	push!(win, g)

	dataPC, gridPC = mkControlGrid([
		("ϵ", Float64, .0, 100.0, 1e-2),
		("κ", Float64, .0, 1.0, 1e-4),
		("δ", Float64, .0, 1e2, 1e-2),
		("α", Float64, .0, 360.0, 1e-1),
		("inv", Bool)
	])
	g[1,1] = gridPC

	buttonStep = @Button("Step")
	g[1,2] = buttonStep
	signal_connect(w -> @async(handlerButtonStep(dataPC)), buttonStep, "clicked")

	g[2,1] = sw = @ScrolledWindow()
	logView  = @TextView()
	push!(sw, logView)
	global PC_logWidgets = (sw, logView)
	setproperty!(logView, :justification, 0)
	setproperty!(logView, :hexpand, true)
	setproperty!(logView, :editable, false)

	showall(win)
end

global PC_logWidgets = false
function showLog(msg)
	if (PC_logWidgets !== false)
		setproperty!(PC_logWidgets[2], :editable, true)
		insert!(PC_logWidgets[2], "\n$msg")
		setproperty!(PC_logWidgets[2], :editable, false)
		adj = getproperty(PC_logWidgets[1], :vadjustment, Gtk.GtkAdjustment)
		setproperty!(adj, :value, getproperty(adj, :upper, Int))
	end
end



function handlerButtonStep(D)
	# choose a solution
	V = Proj[end][end]

	# advance it one step
	W = continuationStep(H, J, V, D["inv"], D["ϵ"], D["κ"], D["δ"], D["α"])

	# write back
	push!(Proj[end], W)
end

# untested
global h = 1.0 #fix
function continuationStep(H, J, V, dir::Bool, ϵ, κ, δ, α)
	global h

	# tangent-vector for matrix A
	function tang(A)
		t = nullspace(A)[:,1]
		copysign(t, det([A; t']))
	end

	# newton corrector step
	function Δ(v)
		Hv, Jv = H(v), J(v)
		Hv, Jv, [Jv; tang(Jv)'] \ [H(v); 0]
	end

	tangJV = tang(J(V))
	v₀ = V + (dir?1:-1) * h * tangJV

	# step size control
	Hv₀, Jv₀, Δv₀ = Δ(v₀)
	v₁ = v₀ - Δv₀
	Hv₁, Jv₁, Δv₁ = Δ(v₁)
	κ′ = norm(Δv₁) / norm(Δv₀)
	δ′ = norm(Δv₀)
	α′ = acos( dot(tang(Jv₁), tang(Jv₀)) )

	f₀ = max(sqrt(κ′ / κ), sqrt(δ′ / δ), α′ / α)
	f = clamp(f₀, .5, 2.0)
	h /= f


	if (f >= 2.0)
		showLog("... h=$h")
		return continuationStep(H, J, V, dir, ϵ, κ, δ, α) #TODO recursive?
	end

	while norm(Hv₁) > ɛ
		Hv₁, Jv₁, Δv₁ = Δ(v₁)
		v₁ = v₁ - Δv₁
	end

	showLog("done!")
	return v₁
end
