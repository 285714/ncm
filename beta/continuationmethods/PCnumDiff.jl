type PCnumDiff <: ContinuationMethod
end


function C(::PCnumDiff)
	# windowPC = @Window("Continuation Settings")

	gridControls = @Grid()
	setproperty!(gridControls, :column_homogeneous, true)
	setproperty!(gridControls, :column_spacing, 10)
	# push!(windowPC, gridControls)

	buttonTest1, buttonTest2 = @Button("Test1"), @Button("Test2")
	scaleTest3 = @Scale(false, 0:100)
	gridControls[1,1] = buttonTest1
	gridControls[2,1] = buttonTest2
	gridControls[1:2,2] = scaleTest3

	# showall(windowPC)

	gridControls
end
