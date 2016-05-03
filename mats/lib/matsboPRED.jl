module matsboPRED
export predCount

# general iterface
# pred must have init::Bool=false keyword argument
# must accept arb number of inputs

# means of putting preds together...
# somehow pipe through v...

function predCount(max::Integer)
	local cnt = 0
	local maxIntern = max
	return function pred(v...; init::Bool=false)
		cnt += 1
		return cnt <= max
	end
end

end
