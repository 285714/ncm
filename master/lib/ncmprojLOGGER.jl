using Gtk.ShortNames

type Logger
	w::Gtk.ScrolledWindowLeaf
	t::Gtk.TextViewLeaf
end

function Logger()
	w = @ScrolledWindow()
	t  = @TextView()
	push!(w, t)
	setproperty!(t, :justification, 0)
	setproperty!(t, :expand, true)
	setproperty!(t, :editable, false)
	return Logger(w,t)
end

import Base.push!
push!(W::Gtk.GtkContainer, L::Logger) = push!(W,L.w)

import Base.write
function write(L::Logger, s)
	setproperty!(L.t, :editable, true)
	insert!(L.t, "$(string(s))\n")
	setproperty!(L.t, :editable, false)
	adj = getproperty(L.w, :vadjustment, Gtk.GtkAdjustment)
	setproperty!(adj, :value, getproperty(adj, :upper, Int))
	return Void
end
