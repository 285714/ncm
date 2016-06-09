using Immerse
using Gtk.ShortNames

#TODO multithreading
#TODO externalize margins, distances etc.

# GUI invariant: Main Window (bifurcation view) is always visible and contains menus etc.

# GUI wrapper
# returns when program is closed
function startGUI()
	windowMain = @Window("Main View", 800, 600, false, true)

	boxImmerse, toolbarImmerse, canvasImmerse = Immerse.createPlotGuiComponents()

	# Main division
	gridMain = @Grid()
	push!(windowMain, gridMain)
	setproperty!(gridMain, :row_spacing, 5)

	# Menu
	menubarMain = @MenuBar()
	gridMain[1:3,1] = menubarMain

	#	File
	menuitemFile = @MenuItem("_File")
	push!(menubarMain, menuitemFile)
	menuFile = @Menu(menuitemFile)

	menuitemFileOpen = @MenuItem("Open")
	push!(menuFile, menuitemFileOpen)

	menuitemFileQuit = @MenuItem("Quit")
	push!(menuFile, menuitemFileQuit)

	#	Cont
	#TODO populate dynamically
	menuitemCont = @MenuItem("_Continuation")
	push!(menubarMain, menuitemCont)
	menuCont = @Menu(menuitemCont)

	menuitemContPC = @MenuItem("PC")
	push!(menuCont, menuitemContPC)

	#	View
	menuitemView = @MenuItem("_View")
	push!(menubarMain, menuitemView)
	menuView = @Menu(menuitemView)

	menuitemViewSSV = @MenuItem("Single Solution View")
	push!(menuView, menuitemViewSSV)

	# Main Plotting Area
	# canvasPlot = @Canvas()
	# setproperty!(canvasPlot, :expand, true)
	# gridMain[1,2] = canvasPlot
	setproperty!(boxImmerse, :expand, true)
	gridMain[1,2] = boxImmerse

	# Control Area
	gridControls = @Grid()
	gridMain[1,3] = gridControls

	setproperty!(gridControls, :column_homogeneous, true)
	setproperty!(gridControls, :column_spacing, 10)
	buttonTest1, buttonTest2 = @Button("Test1"), @Button("Test2")
	scaleTest3 = @Scale(false, 0:100)
	gridControls[1,1] = buttonTest1
	gridControls[2,1] = buttonTest2
	gridControls[1:2,2] = scaleTest3

	showall(windowMain)

	figure(figure(canvasImmerse))
	display(Immerse._display, plot(x=rand(10), y=rand(10)))

	while visible(windowMain) sleep(1/25) end

	return 0
end
