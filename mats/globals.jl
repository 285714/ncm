begin
	# general
	global System, Dimensions
	System, Dimensions = roessler, 3

	# derivative
	global Epsilon
	Epsilon = 1e-2

	# cycle finding
	global TransientIterations, TransientStepSize, SteadyStateIterations, SteadyStateStepSize
	TransientIterations, TransientStepSize = 5000, .1
	SteadyStateIterations, SteadyStateStepSize = 25000, .01

	# initial cycle
	global Samples, InterpPrec
	Samples, InterpPrec = 128, 3
end
