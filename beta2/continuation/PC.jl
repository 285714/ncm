function continuationExec()
	continuationGUI()
end

function continuationGUI()
	win = @Window("ContinuationGUI", 500, 100, false)

	g = @Grid()
	push!(win, g)
	setproperty!(g, :column_spacing, 10)
	setproperty!(g, :row_spacing, 5)
	setproperty!(g, :margin, 10)

	g[1:2,1] = header = @Label("Settings")

	g[1,2] = @Label("ɛ")
	g[2,2] = @SpinButton(0, 1000, 0.00001)

	g[1,3] = @Label("κ")
	g[2,3] = @SpinButton(0, 1000, 0.00001)

	g[1,4] = @Label("δ")
	g[2,4] = @SpinButton(0, 1000, 0.00001)

	g[1,5] = @Label("α")
	g[2,5] = @SpinButton(0, 1000, 0.00001)

	g[1,6] = @Label("Inv")
	g[2,6] = @CheckButton

	g[3,1:6] = sw = @ScrolledWindow()
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



#=
# untested
function continuationStep(pc :: PC, h :: Homotopy, v :: Vector{Float64})
	# tangent-vector for matrix A
	function tang(A)
		t = nullspace(A)
		t * sign(det( [A; t'] ))
	end

	# newton corrector step
	function Δ(v, Jv = J(h)(v))
		[Jv; tang(Jv)'] \ [H(h)(v); 0]
	end

	Jv = J(h)(v)
	v₀ = v + pc.h * pc.direction * tang(Jv)

	# step size control
	Jv₀ = J(h)(v₀)
	Δv₀ = Δ(v₀, Jv₀)
	v₁ = v₀ - Δv₀
	κ′ = norm(Δ(v₁)) / norm(Δv₀)
	δ′ = norm(Δv₀)
	α′ = acos( tang(Jv)' * tang(Jv₀) )[1]

	while norm(H(h)(v₁)) > pc.ɛ
		v₁ = v₁ - Δ(V₁)
	end

	f₀ = max(sqrt(κ′ / pc.κ), sqrt(δ′ / pc.δ), α′ / pc.α)
	f = max(min(f₀, 2), 1/2)
	pc.h = pc.h / f

	if (f < 2)
		log(pc, "step done")
		return v₁
	else
		showLog(pc, "..")
		continuationStep(pc, h, v)
	end
end
=#
