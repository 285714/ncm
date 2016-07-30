module mbObserve
export Observer, observe, unobserve, notify

type Observer
	D::Dict{Any, Vector{Function}}
	Observer() = new( Dict{Any, Vector{Function}}() )
end

function observe(cb::Function, O::Observer, s)
	!haskey(O.D, s) && (O.D[s] = Function[])
	push!(O.D[s], cb)
	return Void
end

function unobserve(cb::Function, O::Observer, s)
	try
		delete!(O.D[s], cb)
		isempty(O.D[s]) && delete!(O.D, s)
	end
	return Void
end

#NOTE not exposed... conceptually dangerous... responsibility for callbacks from different
#	clients can only lie at observer. however clients should know their state, which is
#	not given if Observer can remove their callbacks.
function unobserveall(O::Observer, s)
	O.D[s] = Function[]
	return Void
end

#TODO DRYS
#TODO try... used this way for debugging

function notify(O::Observer, s, v...)
	for cb in get(O.D, s, [])
		# cb(v...)
		try cb(v...) catch e; println(e, "\nSymbol: $s.") end
	end
	return Void
end

function notify(O::Observer, s)
	for cb in get(O.D, s, [])
		# cb()
		try cb() catch e; println(e, "\nSymbol: $s.") end
	end
	return Void
end


end
