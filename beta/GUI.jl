using Gtk.ShortNames

#TODO multithreading
#TODO externalize margins, distances etc.

# GUI invariant: Main Window (bifurcation view) is always visible and contains menus etc.

# GUI wrapper
# returns when program is closed
function startGUI()
	windowMain = @Window("Main View")
	gridMain = @Grid()

	gridControls = @Grid()
	setproperty!(gridControls, :column_homogeneous, true)
	setproperty!(gridControls, :column_spacing, 10)
	buttonTest1, buttonTest2 = @Button("Test1"), @Button("Test2")
	scaleTest3 = @Scale(false, 0:100)

	gridControls[1,1] = buttonTest1
	gridControls[2,1] = buttonTest2
	gridControls[1:2,2] = scaleTest3
	gridMain[4,1] = gridControls
	push!(windowMain, gridMain)

	showall(windowMain)

	while visible(windowMain) sleep(1/25) end
	return 0
end
