
type PC <: ContinuationMethod
	h :: Float64
	ɛ :: Float64
	κ :: Float64
	δ :: Float64
	α :: Float64
	direction :: Int
end

PC() = PC(0.01, 0.01, 1.8, 0.1, 0.51, 1)


# untested
function step(pc :: PC, h :: Homotopy, v :: Vector{Float64})

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
		step(pc, h, v)
	end
end


function showLog(pc :: PC, msg)
	if (PC_logWidgets !== false)
		setproperty!(PC_logWidgets[2], :editable, true)
		insert!(PC_logWidgets[2], "\n$msg")
		setproperty!(PC_logWidgets[2], :editable, false)
		adj = getproperty(PC_logWidgets[1], :vadjustment, Gtk.GtkAdjustment)
		setproperty!(adj, :value, getproperty(adj, :upper, Int))
	end
end
PC_logWidgets = false


macro mkWidgetId(property, args...)
	return :@mkWidget($property, identity, identity, $args)
end

macro mkWidget(property, to, from, args)
	return quote
		let
			t = typeof($from($property))
			(w,p,s) = t == Float64 ? (@SpinButton($args[1], $args[2], $args[3]), :value, :event) :
			          t == Bool    ? (@CheckButton($args[1]), :active, :toggled) : error("no widget for type $t")

			setproperty!(w, p, $from($property))
			signal_connect(w, s) do obj, args...
				v = getproperty(w, p, t)
				$property = $to(v)
			end
			
			w
		end
	end
end

function C(pc :: PC)
	g = @Grid()
	setproperty!(g, :column_spacing, 10)
	setproperty!(g, :row_spacing, 5)
	setproperty!(g, :margin, 10)

	g[1:2,1] = header = @Label("Settings")

	g[1,2] = @Label("ɛ")
	g[2,2] = @mkWidgetId(pc.ɛ, 0, 1000, 0.00001)

	g[1,3] = @Label("κ")
	g[2,3] = @mkWidgetId(pc.κ, 0, 1000, 0.00001)

	g[1,4] = @Label("δ")
	g[2,4] = @mkWidgetId(pc.δ, 0, 1000, 0.00001)

	g[1,5] = @Label("α")
	g[2,5] = @mkWidgetId(pc.α, 0, 1000, 0.00001)

	g[1:2,6] = @mkWidget(pc.direction, x -> x ? -1 : 1, x -> x < 0, ["inverted direction"])

	g[3,1:6] = sw = @ScrolledWindow()
	logView  = @TextView()
	push!(sw, logView)
	global PC_logWidgets = (sw, logView)
	setproperty!(logView, :justification, 0)
	setproperty!(logView, :hexpand, true)
	setproperty!(logView, :editable, false)

	g
end

