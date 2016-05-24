global ax1, ax2, distance

function initPlot3D()
		global ax1, ax2, distance
		distance = Float64[]
    	ion(); figure()
		global ax1 = subplot2grid((4,1), (0, 0), rowspan=3, projection="3d"); plot([.0],[.0],[.0])
		global ax2 = subplot2grid((4,1), (3, 0))
end

function callbackPlot3D(V, H, J)
	global ax1, ax2, distance
	local y = toTraj(V[1:end-2])

	sca(ax1); ax1["lines"][end]["set_color"]("gray"); hold(true); plot(y[:,1], y[:,2], y[:,3], color="r")

	push!(distance, sumabs(H))
	sca(ax2); hold(false); plot(distance)

	draw(); waitforbuttonpress()
end

function plotLast(V)
	global ax1
	local C,ω,N
	C,ω = rToC(V)
	N = 2*(size(C,1)-1)
	sca(ax1); hold(true); plot(irfft(C[:,1],N), irfft(C[:,2],N), irfft(C[:,3],N), color="g")
end






function initPlotResidual()
		global ax1, ax2, distance
		distance = Float64[]
    ion(); figure()
		global ax1 = subplot2grid((4,1), (0, 0), rowspan=3); plot([.0],[.0])
		global ax2 = subplot2grid((4,1), (3, 0))
end

function callbackPlotResidual(V, H, J)
	global ax1, ax2, distance

  sca(ax1); ax1["lines"][end]["set_color"]("gray"); hold(true);	plot(H, color="r")

	push!(distance, sumabs(H))
	sca(ax2); hold(false); plot(distance)

	draw()
	waitforbuttonpress()
end




initIncLoad() = begin ion(); subplot(111,projection="3d"); show(); sleep(.01) end
callbackIncLoad(::Any,::Any,w) = if w!=Void
	plot(w[:,1], w[:,2], w[:,3])
	draw(); sleep(.01)
end





function initPC()
	global pl
	for p in pl
		p["set_xdata"](reverse(p["get_xdata"]()))
		p["set_ydata"](reverse(p["get_ydata"]()))
	end

	hold(true)
	ion()
end

# function callbackPC(V)
# 	global pl
# 	local v,ℵ,A
#
# 	v,ℵ = toTraj(V[1:end-2]), V[end]
# 	A = Float64[]
#
# 	f(a,b) = -a[1]/(-a[1]+b[1])*a + b[1]/(-a[1]+b[1])*b
#
# 	for i in 1:size(v,1)-1
# 		v[i,1] ≤ 0 && v[i+1,1] > 0 && push!(A,norm(f(v[i,:],v[i+1,:])))
# 	end
# 	v[end,1] ≤ 0 && v[1,1] > 0 && push!(A,norm(f(v[end,:],v[1,:])))
#
# 	while length(pl) < length(A); push!(pl, plot([],[])[1]) end
# 	for i in 1:length(A)
# 		local x,y; x,y = pl[i]["get_xdata"](), pl[i]["get_ydata"]()
# 		push!(x,ℵ); push!(y,A[i])
# 		pl[i]["set_xdata"](x); pl[i]["set_ydata"](y)
# 	end
# 	gca()["relim"](); gca()["autoscale"]()
# 	sleep(.01)
#
# 	return Void
# end




function callbackPC(V)
	global pl, InternSamples
	local P,ℵ,A

	P,ℵ = toTrajInterp(V[1:end-2]), V[end]
	A = Float64[]

	local t = linspace(.0,2pi,InternSamples)
	for (a,b) in zip(t, t+t[2])
		P[1](a) ≤ 0 ≤ P[1](b) || continue
		local c = bisection(P[1],a,b)
		push!(A, norm(map(f->f(c), P[2:end])))
	end

	while length(pl) < length(A); push!(pl, plot([],[])[1]) end
	for i in 1:length(A)
		local x,y; x,y = pl[i]["get_xdata"](), pl[i]["get_ydata"]()
		push!(x,ℵ); push!(y,A[i])
		pl[i]["set_xdata"](x); pl[i]["set_ydata"](y)
	end
	gca()["relim"](); gca()["autoscale"]()
	sleep(.01)

	return Void
end
