module mbObserve
export Observer, observe, unobserve, notify

type Observer
	D::Dict{Symbol, Vector{Function}}
	Observer() = new( Dict{Symbol, Vector{Function}}() )
end

function observe(cb::Function, O::Observer, s::Symbol)
	!haskey(O.D, s) && (O.D[s] = Function[])
	push!(O.D[s], cb)
	return Void
end

function unobserve(cb::Function, O::Observer, s::Symbol)
	try delete!(O.D[s], cb) end
	return Void
end

#NOTE not exposed... conceptually dangerous... responsibility for callbacks from different
#	clients can only lie at observer. however clients should know their state, which is
#	not given if Observer can remove their callbacks.
function unobserveall(O::Observer, s::Symbol)
	O.D[s] = Function[]
	return Void
end

function notify(O::Observer, s::Symbol, v...)
	for cb in get(O.D, s, [])
		cb(v...)
	end
	return Void
end

function notify(O::Observer, s::Symbol)
	for (s,F) in O.D
		for cb in F
			cb()
		end
	end
	return Void
end


end
