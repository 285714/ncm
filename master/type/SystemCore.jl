# responsible for the interpretation of solutions
abstract SystemCore

# usually displays a GUI
Base.show(::SystemCore) = error("implement!")

# return  R^(N+1)->R^N  /  R^(N+1)->R^(Nx(N+1))  functions
H(::SystemCore) = error("implement!")
J(::SystemCore) = error("implement!")

#TODO replace observer with required functions in system core?
#	interface should be clear...
