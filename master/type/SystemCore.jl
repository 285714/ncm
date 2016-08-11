# responsible for the interpretation of solutions
"""
An abstract type that requires to implement the basic functions needed for path following:
Homotopy `H`, its jacobian `J` (and `show` to display a GUI)
"""
abstract SystemCore

# usually displays a GUI
Base.show(::SystemCore) = error("implement!")

# return  R^(N+1)->R^N  /  R^(N+1)->R^(Nx(N+1))  functions
H(::SystemCore) = error("implement!")
J(::SystemCore) = error("implement!")

#TODO replace observer with required functions in system core?
#	interface should be clear...
