#TODO multithreading
#TODO externalize margins, distances etc.

# GUI invariant: Main Window (bifurcation view) is always visible and contains menus etc.

# GUI wrapper
# returns when program is closed
function startGUI()
	global windowMain = @Window("Main View", 800, 600, false, true)

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

	menuitemFileNew = @MenuItem("New")
	push!(menuFile, menuitemFileNew)

	menuitemFileOpen = @MenuItem("Open")
	push!(menuFile, menuitemFileOpen)

	menuitemFileQuit = @MenuItem("Quit")
	push!(menuFile, menuitemFileQuit)

	#=
	#	Cont
	#TODO populate dynamically
	menuitemCont = @MenuItem("_Continuation")
	push!(menubarMain, menuitemCont)
	menuCont = @Menu(menuitemCont)

	menuitemContPC = @MenuItem("PC")
	push!(menuCont, menuitemContPC)
	=#

	#	View
	menuitemView = @MenuItem("_View")
	push!(menubarMain, menuitemView)
	menuView = @Menu(menuitemView)

	menuitemViewBV = @MenuItem("Bifurcation View")
	push!(menuView, menuitemViewBV)

	menuitemViewSSV = @MenuItem("Single Solution View")
	push!(menuView, menuitemViewSSV)

	# Main Plotting Area
	# canvasPlot = @Canvas()
	# setproperty!(canvasPlot, :expand, true)
	# gridMain[1,2] = canvasPlot
	setproperty!(boxImmerse, :expand, true)
	gridMain[1,2] = boxImmerse

	# Control Area ~> now in separate window
	# depending on continuation method
	#gridMain[1,3] = C(Proj.Cont)

	# Signals
	signal_connect(menuitemFileNewActivateHandler, menuitemFileNew, "activate")
	signal_connect(menuitemFileOpenActivateHandler, menuitemFileOpen, "activate")
	signal_connect(menuitemViewBVActivateHandler, menuitemViewBV, "activate")
	signal_connect(menuitemViewSSVActivateHandler, menuitemViewSSV, "activate")

	showall(windowMain)

	figure(canvasImmerse)
	display(plot(y=rand(10), Geom.line()))
	return 0
end


#TODO disable unavailable items

# Handlers
#TODO singleton handling
menuitemViewSSVActivateHandler(widget) = C(Proj.Hom)
menuitemViewBVActivateHandler(widget) = C(Proj.Cont)

#TODO unsaved bla
function menuitemFileNewActivateHandler(widget)
	windowNew = @Window("New", 200, 100, false)

	gridNew = @Grid()
	push!(windowNew, gridNew)

	comboCont = @ComboBoxText()
	for c in subtypes(ContinuationMethod)  push!(comboCont, string(c)) end
	gridNew[1,1] = comboCont

	comboHom = @ComboBoxText()
	for c in subtypes(Homotopy)  push!(comboHom, string(c)) end
	gridNew[1,2] = comboHom

	showall(windowNew)
end



function menuitemFileOpenActivateHandler(widget)
	include(open_dialog("Choose a system file.", windowMain))
end
