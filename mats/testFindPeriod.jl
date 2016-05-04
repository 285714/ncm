using PyPlot

include("lib/matsboRK.jl")
include("lib/matsboUTIL.jl")
include("lib/matsboPRED.jl")

include("lib/ncmprojFINDINITIALDATA.jl")

include("systems.jl")

#careful: lorenz diverges for too large TransientStepSizes...
data, cyc = findRepresentativeCycle(roessler, .0, rand(3), 5000, .1, 25000, .01)
data, cyc = findRepresentativeCycle(lorenz, .0, rand(3), 5000, .01, 25000, .01)

figure(); subplot(111,projection="3d")
plot(data[:,1], data[:,2], data[:,3], color="b")
plot(cyc[:,1], cyc[:,2], cyc[:,3], color="r")

# periodic regression? - optimize fourier coefficients to minimize dist to points
# generally forget sequence of points?
# thin out trajectories for AC?
