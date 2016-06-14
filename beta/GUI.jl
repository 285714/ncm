#TODO multithreading
#TODO externalize margins, distances etc.

# GUI invariant: Main Window (bifurcation view) is always visible and contains menus etc.

# GUI wrapper
# returns when program is closed


macro mkWidgetId(property, args...)
	return :@mkWidget($property, identity, identity, $args)
end

macro mkWidget(property, to, from, args)
	return quote
		let
			t = typeof($from($property))
			(w,p,s) = t <: Float64 ? (@SpinButton($args[1], $args[2], $args[3]), :value, :event) :
			          t == Bool   ? (@CheckButton($args[1]), :active, :toggled) : error("no widget for type $t")

			setproperty!(w, p, $from($property))
			signal_connect(w, s) do obj, args...
				v = getproperty(w, p, t)
				$property = $to(v)
			end
			
			w
		end
	end
end



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
	isdefined(:Proj) && (gridMain[1,3] = continuationControlArea() |> showall) # has to be triggered again

	# Signals
	signal_connect(menuitemFileNewActivateHandler, menuitemFileNew, "activate")
	signal_connect(menuitemFileOpenActivateHandler, menuitemFileOpen, "activate")
	signal_connect(menuitemViewBVActivateHandler, menuitemViewBV, "activate")
	signal_connect(menuitemViewSSVActivateHandler, menuitemViewSSV, "activate")

	showall(windowMain)

	figure(canvasImmerse)
	display(plot(x=rand(100), y=rand(100)))

	return 0
end


#TODO disable unavailable items

# Handlers
#TODO singleton handling
menuitemViewSSVActivateHandler(widget) = C(Proj.Hom)

let
	windowCont = false
	oldMethod = false
	global menuitemViewBVActivateHandler = function(widget)
		if windowCont !== false
			if oldMethod != typeof(Proj.Cont)
				empty!(windowCont)
				push!(windowCont, C(Proj.Cont)) |> showall
				oldMethod = typeof(Proj.Cont)
			end
			present(windowCont)
			return
		end

		windowCont = @Window("Continuation Settings", 500, 220, false) |> C(Proj.Cont) |> showall
		signal_connect(_ -> windowCont = false, windowCont, :destroy)
		oldMethod = typeof(Proj.Cont)
	end


	controlSettings = [1.0] # [steps]
	usedContMethods = Dict{AbstractString,ContinuationMethod}()
	global continuationControlArea = function()
		b = @Box(:h)

		setproperty!(b, :spacing, 20)
		setproperty!(b, :margin, 5)

		b |>
			(tools = @Toolbar()) |>
			@SeparatorToolItem() |>

			@Label("Number of steps") |>
			@mkWidgetId(controlSettings[1], 1, 10000, 1) |>
			@SeparatorToolItem() |>

			@Label("Method") |>
			(comboCont = @ComboBoxText())

		run = @ToolButton("Run")
		setproperty!(run, :label, "Run")
		setproperty!(run, :icon_name, "media-playback-start")
		insert!(tools, 0, run)
		signal_connect(run, :clicked) do obj, args...
			n = controlSettings[1]
			runContinuationMethod(Proj.Cont, Proj.Hom, [0.0], (c,h,v,n′) -> n′ < n)
		end

		pause= @ToolButton("Pause")
		setproperty!(pause, :label, "Pause")
		setproperty!(pause, :icon_name, "media-playback-pause")
		insert!(tools, 0, pause)
		signal_connect(pause, :clicked) do obj, args...
			println("stopping..")
		end

		contMethods = subtypes(ContinuationMethod)
		for method in contMethods  push!(comboCont, string(method)) end
		signal_connect(comboCont, :changed) do obj, args...
			method = contMethods[getproperty(comboCont, :active, Int) + 1]
			if haskey(usedContMethods, string(method))
				methodInst = get(usedContMethods, string(method), 0)
			else
				methodInst = try method() catch error("no trivial constructor for $method") end
				push!(usedContMethods, string(method), methodInst)
			end
			Proj.Cont = methodInst
			windowCont === false || menuitemViewBVActivateHandler(false)
		end
		push!(usedContMethods, string(typeof(Proj.Cont)), Proj.Cont)
		setproperty!(comboCont, :active, find(x -> x == typeof(Proj.Cont), contMethods)[1] - 1)

		b
	end
end


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
