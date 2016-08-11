#TODO convert to real MVC, unify interface, plugin style, fix grid interface
#TODO additional abstraction layer encapsulating an element plus data and handler

"""
    ctrl(D, x)
Dictionary `D`, Tuple `x=(name::String, ::DataType, init, v...)`; used by mkControlGrid
"""
function ctrl(D::Dict{AbstractString, Any}, x)
	local tmp = Void
	if x[2] == Int
		tmp = @SpinButton(x[4][1], x[4][end], step(x[4]))
		setproperty!(tmp, :value, x[3])
		signal_connect(tmp, "value_changed") do w
			D[x[1]] = getproperty(w, :value, x[2])
		end
		signal_emit(tmp, "value_changed", Void)
	elseif x[2] == Float64
		tmp = @Scale(:h, x[4][1], x[4][end], step(x[4]))
		local adj = @Adjustment(tmp)
		setproperty!(adj, :value, x[3])
		signal_connect(adj, "value_changed") do w
			D[x[1]] = getproperty(w, :value, x[2])
		end
		signal_emit(tmp, "value_changed", Void)
	elseif x[2] == Bool
		tmp = @CheckButton()
		setproperty!(tmp, :active, x[3])
		signal_connect(tmp, "toggled") do w
			D[x[1]] = getproperty(w, :active, x[2])
		end
		signal_emit(tmp, "toggled", Void)
	end

	if tmp == Void
		warn("No handling for type $(x[2]) defined.")
		return Void
	end

	D[x[1]] = x[3]

	setproperty!(tmp, :hexpand, true)
	setproperty!(tmp, :vexpand, false)
	return tmp
end

function button(handler::Function, label::AbstractString)
	B = @Button(label)
	signal_connect(B, "clicked") do w
		@schedule handler(w)
	end
	return B
end


"""
    mkControlGrid(D, C)
Creates a grid of controls with labels, handlers and encapsulated storage (Dictionary `D`)
`c` in `C` is Tuple `(name::String, ::DataType, init, v...)`
"""
function mkControlGrid(D::Dict{AbstractString, Any}, C)
	Control = map(C) do X
		reduce([], X) do A,a
			isa(a, Gtk.GtkWidget) && return [A; a]
			a == Void && return [A; Void]
			return [A; @Label(a[1]); ctrl(D, a)]
		end
	end

	G = @Grid()
	setproperty!(G, :column_spacing, 10)
	setproperty!(G, :row_spacing, 6)
	setproperty!(G, :margin, 6)
	setproperty!(G, :column_homogeneous, true)
	# setproperty!(G, :row_homogeneous, true)


	S = reduce((A,a)->lcm(A,length(a)), 1, Control)
	for i in 1:length(Control)
		l = S÷length(Control[i])
		for j in 1:length(Control[i])
			Control[i][j] == Void && continue
			G[(j-1)*l + (1:l), i] = Control[i][j]
		end
	end

	return G
end
