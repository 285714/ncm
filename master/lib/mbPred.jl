__precompile__()

module mbPred
export predCount, predEps, ∧

# general iterface:
#		pred are 'until'-predicates, becoming false when a loop should stop
#		functions here are constructors for predicates
#		pred must have init::Bool=false keyword argument
#			on init==true initialize and return true, afterwards behaviour should be same as when newly created
#		? must accept arb number of inputs

# means of putting preds together...
# somehow pipe through v...

# and-connect predicates
∧(p...) = function pred(v...; init::Bool=false)
	return foldl( (A,a) -> A&&a, map( f->f(v..., init=init), p ) )
end

function predCount(max::Integer)
	local cnt = 0
	return function pred(v...; init::Bool=false)
		if init
			cnt = 0
			return true
		end
		cnt += 1
		return cnt <= max
	end
end


function predEps(ϵ,i=1)
	return function pred(v...; init::Bool=false)
		if init return true end
		return norm(v[i]) ≥ ϵ
	end
end

end
