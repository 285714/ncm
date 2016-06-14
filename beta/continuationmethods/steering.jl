
# use Tasks..
function runContinuationMethod(cm :: ContinuationMethod, h :: Homotopy, v :: Vector{Float64}, p :: Function)
	n = 0
	while p(cm, h, v, n) !== false
		println("step..")
		# v = step(cm, h, v)
		n += 1
	end
end




