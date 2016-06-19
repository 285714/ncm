using matsboNWTN
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
		("α", Float64, .0, 180.0, 1e-1),
		("inv", Bool)
	])
	g[1,1] = gridPC

	buttonStep = @Button("Step")
	g[1,2] = buttonStep
	signal_connect(w -> handlerButtonStep(dataPC), buttonStep, "clicked") #TODO async

	buttonRun = @ToggleButton("Run")
	g[1,3] = buttonRun
	signal_connect(w -> handlerButtonRun(w, dataPC), buttonRun, "toggled")
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

function handlerButtonRun(w, D)
	@async while getproperty(w, :active, Bool)
			handlerButtonStep(D)
			yield()
	end
end

#TODO dir...
function handlerButtonStep(D)
	# choose a solution
	V = Proj[end][end]

	# advance it one step
	W = continuationStep(H, J, V, D["inv"], D["ϵ"], D["κ"], D["δ"], D["α"])

	# write back
	push!(Proj[end], W)

	plotBifurcation()
end

# untested
global h = 1.0 #fix
function continuationStep(H, J, V, dir::Bool, ϵ, κ, δ, α)
	# tangent-vector for matrix A
	function tang(A)
		t = A \ [zeros(size(A,1)-1); 1]
		copysign(t / norm(t), det([A; t']))
	end

	# newton corrector step
	function Δ(v)
		Hv, Jv = H(v), J(v)
		Hv, Jv, Jv \ Hv
	end

	local v₁, Hv₁
	while true
		# showLog("... h=$h") #TODO fix... blocks execution at that point
		v₀ = V + (dir?h:-h) * tang(J(V))

		# step size control
		Hv₀, Jv₀, Δv₀ = Δ(v₀)
		v₁ = v₀ - Δv₀
		Hv₁, Jv₁, Δv₁ = Δ(v₁)
		v₁ = v₁ - Δv₁
		κ′ = norm(Δv₁) / norm(Δv₀)
		δ′ = norm(Δv₀)
		α′ = acos( dot(tang(Jv₁), tang(Jv₀)) )

		f = clamp(max(sqrt(κ′/κ), sqrt(δ′/δ), α′/α), .5, 2.0)
		global h /= f

		f < 2.0 && break
	end

	println("$h")
	newton(H, J, v₁, predEps(ϵ))

	#showLog("done!")
	return v₁
end
