include("$(pwd())/lib/ncmprojFINDINITIALDATA.jl")
include("$(pwd())/lib/ncmprojMKCONTROLGRID.jl")

type Galerkin <: SystemCore
	parent::Session
	f::Function
	projection::Function
	H::Function #derived
	J::Function #derived
	dataGUI::Dict{AbstractString, Any}
	Galerkin(parent,f,f′,projection) = new(parent,f,projection,fToH(f),f′ToJ(f′),Dict{AbstractString, Any}())
	Galerkin(parent,f,H,J,projection) = new(parent,f,projection,H,J,Dict{AbstractString, Any}())
end

#TODO dimensions

H(G::Galerkin) = G.H
J(G::Galerkin) = G.J

function Base.show(G::Galerkin)
	windowGalerkin = @Window("Galerkin Controls", 256, 100, false, true)
	setproperty!(windowGalerkin, :resizable, false)
	# setproperty!(windowGalerkin, :type_hint, Gtk.GdkWindowTypeHint.TOOLBAR)

	gridGalerkin = @Grid()
	setproperty!(gridGalerkin, :column_spacing, 5)
	setproperty!(gridGalerkin, :row_spacing, 5)
	push!(windowGalerkin, gridGalerkin)

	buttonFindInitial = @ToggleButton("Find Initial")
	local taskInitial
	signal_connect(buttonFindInitial, "toggled") do w
		if getproperty(w, :active, Bool)
			taskInitial = @schedule try handlerFindInitialData(w, G)
				finally @schedule setproperty!(w, :active, false) end
		else
			isa(taskInitial, Task) && !istaskdone(taskInitial) && schedule(taskInitial, InterruptException(); error=true)
		end
	end

	gridGalerkin[1:2,1] = mkControlGrid(G.dataGUI, Array[
		[ ("Trans. Iterations", Int, 4000, 1000:1000:1e8), Void, Void ],
		[ ("Trans. StepSize", Float64, .01, 1e-3:.01:1.0), ("SS StepSize", Float64, .01, 1e-3:.01:1.0) ],
		[ ("Max. period", Int, 30, 1:128), ("Intersections", Int, 120, 1:128) ],
		[ ("c₀", Float64, 350, -500:.01:500), ("m", Int, 64, 8:4096) ],
		[ buttonFindInitial ],
		[ ("Samples", Int, 64, 8:4096), ("Factor", Float64, 1.0, .0:.1:8.0) ],
		[ button(w->handlerResample(w, G), "Process") ]
	])

	showall(windowGalerkin)
end

function handlerResample(w, G)
	setproperty!(w, :sensitive, false)
	try
		#handle !isa Solution
		C = resample(G.parent.P.activeSolution, G.dataGUI["Samples"], G.dataGUI["Factor"])
		ω,ℵ = C[end-1:end]

		Htmp(V) = G.H([V; ℵ])
		Jtmp(V) = G.J([V; ℵ])[:, 1:end-1]
		tmp = newton(Htmp, Jtmp, C[1:end-1], predCount(10) ∧ predEps(1e-10))

		setActiveSolution(G.parent.P, [tmp; ℵ])
	finally
		setproperty!(w, :sensitive, true)
	end
	return Void
end


function handlerFindInitialData(w, G)
	try
		#TODO function/macro bringIntoScope(D::Dict)
		# TIters, SSIters, TStepSize, SSStepSize, Periods, m, c₀ = map(x->G.dataGUI[x], ["Trans. Iterations", "SS Iterations", "Trans. StepSize", "SS StepSize", "Periods", "m", "c₀"])
		TIters, TStepSize, SSStepSize, maxCycles, nIntersections, m, c₀ = map(x->G.dataGUI[x], ["Trans. Iterations", "Trans. StepSize", "SS StepSize", "Max. period", "Intersections", "m", "c₀"])

		# dataT, dataSS, P = findCycle((t,v)->f(t,[v;c₀]), .0, rand(3), TIters, TStepSize, SSIters, SSStepSize)
		# cyc,ω = prepareCycle(dataSS, SSStepSize, P; fac=Periods)
		cyc, P, ω = findCyclePoincare((t,v)->G.f(t,[v;c₀]), rand(3),
			nIntersections=nIntersections, maxCycles=maxCycles, sampleSize=m,
			transientIterations=TIters, transientStepSize=TStepSize,
			steadyStateStepSize=SSStepSize)

		C = wrap(rfft(cyc, [1]), ω, c₀)
		C = resample(C, m)
		C = C[1:end-1]

		Htmp(V) = G.H([V; c₀])
		Jtmp(V) = G.J([V; c₀])[:, 1:end-1]
		C = newton(Htmp, Jtmp, C, predCount(10) ∧ predEps(1e-10))

		setActiveSolution(G.parent.P, [C; c₀])
	catch e
		println(e)
	end
	return Void
end

function unwrap(V)
	isa(V, Solution) && (V=V.data)
	ω,ℵ = V[end-1:end]
	m = (length(V)-5)÷6
	W = reshape(V[1:end-2], 2m+1, 3)
	W = complex(W[1:m+1,:], [zeros(1, size(W,2)); W[m+2:end,:]])
	return W,ω,ℵ
end

function wrap(V,ω,ℵ)
	return [vec(vcat(real(V), imag(V[2:end,:]))); ω; ℵ]
end

# array of functions interpolating a component each
function interps(W)
	m = size(W,1)-1
	fs = mapslices(W, [1]) do w
		interpolateTrigonometric(real(w[1])/2m, real(w[2:end])/m, -imag(w[2:end])/m)
	end
	return vec(fs)
end

# complete vectorized interpolated trajectory
#TODO this is trivially parallelizable...
function interp(W)
	fs = interps(W)
	return t -> mapreduce(g -> g(t), hcat, fs)
end

function projection(G::Galerkin, V)
	W,ω,ℵ = unwrap(V)
	m = size(W,1)-1
	return G.projection(interps(W), ℵ, 2m)
end


# takes real F-coefficients  V  , returns resampled version
function resample(V, m::Int, scale=1.0)
	W,ω,ℵ = unwrap(V)
	f = interp(W)
	T = linspace(.0,2pi*scale,2*m+2)[1:end-1]
	wrap(rfft(f(T),[1]), ω/scale, ℵ)
end


function fToH(f)
	return function H(V::Vector{Float64})
		m = length(V-5)÷6
		anchor = sum(V[1:m+1]) #TODO hm...

		ω,ρ = V[end-1:end]
		V = reshape(V[1:end-2],2m+1,3)
		V = complex(V[1:m+1,:], [zeros(3)'; V[m+2:end,:]])

		defect = rfft(mapslices(x->f(0,[x;ρ]), irfft(V, 2m, [1]), [2]), [1]) - im*ω*(0:m) .* V
		defect = vec([ real(defect); imag(defect)[2:end,:] ])

		return [defect; anchor]
	end
end

function f′ToJ(f′)
	return function J(V::Vector{Float64})
		m = length(V-5)÷6
		ω,ρ = V[end-1:end]
		V = reshape(V[1:end-2],2m+1,3)
		V = complex(V[1:m+1,:], [zeros(3)'; V[m+2:end,:]])

		T = irfft(V, 2m, [1])
		M = Array{Float64}(size(T,1), size(T,2), size(T,2)+1) #n, comp, deriv
		for i in 1:size(T,1) M[i,:,:] = f′(0,[T[i,:]'; ρ]) end

		M′ = [ rCCD(rfft(M[:,i,j])) for i in 1:size(T,2), j in 1:size(T,2) ]

		for i in 1:size(T,2)
			D = diagm(ω*(0:m))
			M′[i,i][m+2:end,1:m+1] -= D[2:end,:]
			M′[i,i][1:m+1, m+2:end] += D[:,2:end]
		end

		M′ = reducedim(vcat, M′, [1], Array{Float64}(0,size(M′[1,1],2)))
		M′ = reducedim(hcat, M′, [2], Array{Float64}(size(M′[1,1],1),0))[1]

		ddω = -im*(0:m).*V
		ddω = vec([real(ddω); imag(ddω)[2:end,:]])

		ddρ = reduce(hcat, [ rfft(M[:,i,size(T,2)+1]) for i in 1:size(T,2) ])
		ddρ = vec([real(ddρ); imag(ddρ)[2:end,:]])

		return [
			M′ ddω ddρ
			[ones(1,m+1) zeros(1,5m+2)] 0 0
		]
	end
end




# Jacobian of circular convolution of coefficients of real functions in appropriate format.
# ∂/∂X cconv(X,Y) = ∂/∂X F(x⋅y) = rCCD(Y)
rCCD(C) = rCCD(real(C[1]), real(C[2:end]), imag(C[2:end]))
function rCCD(V₀,Vᵣ,Vᵢ)
	m = length(Vᵣ)
	I1 = [ mod(i-j, 2m+1)+1 for i in 0:m, j in 0:m ]
	I2 = [ mod(i+j, 2m+1)+1 for i in 0:m, j in 0:m ]
	Wᵣ = [V₀; Vᵣ; Vᵣ[end:-1:1]]
	Wᵢ = [.0; Vᵢ; -Vᵢ[end:-1:1]]

	local rtn =  [
		(Wᵣ[I1]+Wᵣ[I2])					(-Wᵢ[I1]+Wᵢ[I2])[:,2:end]
		(Wᵢ[I1]+Wᵢ[I2])[2:end,:]		(Wᵣ[I1]-Wᵣ[I2])[2:end,2:end]
	] / (2m)
	rtn[:,1] /= 2
	return rtn
end
