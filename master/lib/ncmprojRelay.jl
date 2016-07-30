# relay execution of certain functions on type  T  to one of its fields  f  .
# kind of a workaround for missing inheritance
#TODO rethink this... might fuck up interface e.g. in conjuction with observer
macro relay(T′, f′)
	T,f = eval(T′), eval(f′)	# use esc instead?
	@assert f in fieldnames(T)
	U = fieldtype(T, f)

	return quote
		Base.setindex!(x::$T, v...) 	= setindex!(x.$f, v...)
		Base.getindex(x::$T, v...) 		= getindex(x.$f, v...)
		Base.start(x::$T) 				= start(x.$f)
		Base.next(x::$T, state) 		= next(x.$f, state)
		Base.done(x::$T, state) 		= done(x.$f, state)
		Base.convert(::Type{$U}, x::$T)	= x.$f
		Base.push!(x::$T, items...)		= push!(x.$f, items...)
		Base.pop!(x::$T, items...)		= pop!(x.$f, items...)
		Base.unshift!(x::$T, items...)	= unshift!(x.$f, items...)
		Base.shift!(x::$T, items...)	= shift!(x.$f, items...)
		Base.length(x::$T)				= length(x.$f)
		Base.endof(x::$T)				= endof(x.$f)
		return Void
	end
end
