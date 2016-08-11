"""
Continuation method implementations extend  `ContinuationMethod`
and the corresponding `show` and `step`.
Responsible for changing the project itself
"""
abstract ContinuationMethod

# generic show function. usually displays a GUI
Base.show(::ContinuationMethod) = error("implement!")

# step function. advances activeSolution by one step.
# atm called only by the cont gui.
# exposed so it can be called eg from a bifurcation view.
Base.step(::ContinuationMethod) = error("implement!")
