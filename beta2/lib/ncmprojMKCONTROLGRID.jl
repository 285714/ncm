#TODO convert to real MVC, unify interface, plugin style, fix grid interface
#TODO additional abstraction layer encapsulating an element plus data and handler

# creates a grid of controls with labels, handlers and encapsulated storage
# c in C is Tuple (name::String, ::Type, init, v...)
function mkControlGrid(D::Dict{AbstractString, Any}, C, cols=1)
	for x in C; D[x[1]] = x[3] end

	local Control = map(C) do x
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
		tmp == Void && error("No handling for type $(x[2]) defined.")
		setproperty!(tmp, :hexpand, true)
		return tmp
	end

	local G = @Grid()
	setproperty!(G, :column_spacing, 10)
	setproperty!(G, :row_spacing, 6)
	setproperty!(G, :margin, 6)

	I,R = divrem(length(C),cols)
	for j in 1:cols, i in 1:I + (j<=R?1:0)
		local idx = (j-1)*I + clamp(j-1,0,R) + i
		G[2j-1, i] = @Label(string(C[idx][1]))
		G[2j, i] = Control[idx]
	end

	return G
end
