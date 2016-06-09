type PC <: ContinuationMethod
end


function C(::PC)
	windowPC = @Window("Continuation Settings", 300, 100, false)

	gridControls = @Grid()
	setproperty!(gridControls, :column_homogeneous, true)
	setproperty!(gridControls, :column_spacing, 10)
	push!(windowPC, gridControls)

	buttonTest1, buttonTest2 = @Button("Test1"), @Button("Test2")
	scaleTest3 = @Scale(false, 0:100)
	gridControls[1,1] = buttonTest1
	gridControls[2,1] = buttonTest2
	gridControls[1:2,2] = scaleTest3

	showall(windowPC)
end
