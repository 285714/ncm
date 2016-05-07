# general
global System, Dimensions
System, Dimensions = lorenz, 3

# derivative
global Epsilon
Epsilon = 1e-8

# cycle finding
global TransientIterations, TransientStepSize, SteadyStateIterations, SteadyStateStepSize
TransientIterations, TransientStepSize = 5000, .01
SteadyStateIterations, SteadyStateStepSize = 25000, .005

# initial cycle
global Samples
Samples = 512
