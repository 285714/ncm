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

  local C,ω,N
  C,ω = rToC(V)
	N = 2*(size(C,1)-1)

  sca(ax1); ax1["lines"][end]["set_color"]("gray"); hold(true); plot(irfft(C[:,1],N), irfft(C[:,2],N), irfft(C[:,3],N), color="r")

	push!(distance, sumabs(H))
	sca(ax2); hold(false); plot(distance)

	draw(); sleep(.05)
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
