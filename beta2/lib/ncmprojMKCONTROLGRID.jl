# creates a grid of controls with labels, handlers and encapsulated storage
function mkControlGrid(C, cols=1)
	local Data = [ x[1] => zero(x[2]) for x in C ]

	local Control = map(C) do x
		local tmp = Void
		if x[2] == Int
			tmp = @SpinButton(x[3:end]...)
			signal_connect(tmp, "value_changed") do w
				Data[x[1]] = getproperty(w, :value, x[2])
			end
		elseif x[2] == Float64
			tmp = @Scale(:h, x[3:end]...)
			local adj = @Adjustment(tmp)
			signal_connect(adj, "value_changed") do w
				Data[x[1]] = getproperty(w, :value, x[2])
			end
		elseif x[2] == Bool
			tmp = @CheckButton(x[3:end]...)
			signal_connect(tmp, "toggled") do w
				Data[x[1]] = getproperty(w, :value, x[2])
			end
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

	return Data, G
end
