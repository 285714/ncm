
abstract ContinuationMethod

# encapsulates the state of the Predictor-Corrector method
type PredictorCorrector <: ContinuationMethod	
	h :: Real		# Stepsize
	ɛ :: Real		# Newton ɛ
	κ :: Real
	δ :: Real
	α :: Real
	direction :: Int

	v				# current Position
	n :: Int		# number of steps made
end

PredictorCorrector() = PredictorCorrector(0.01, 0.01, 1.8, 0.1, 0.51, 1, (), 0)


type Broyden
	# TODO
end

Broyden() = Broyden()


function initialize end
function step end
function unwrap(x :: Any) x end
function wrap!(x :: Any, V) V end

macro wrapped(obj, f)
	return :(v -> $f(wrap!($obj, v)))
end


function cm(method :: ContinuationMethod,
	        H :: Function,
			DH :: Function,
			callback :: Function,
			obj = ())

	if obj != ()
		initialize(method, H, DH, unwrap(obj))
	end

	while callback(method) !== false
		step(method, H, DH)
	end
end


export PredictorCorrector, Broyden, cm, wrapped



# implementations for simple predictor-corrector method

function initialize(method :: PredictorCorrector, H, DH, v₀)
	# TODO: initial newton
	method.v = v₀
end

function step(method :: PredictorCorrector, H, DH)
	method.n += 1

	function tang(A)
		t = nullspace(A)
		t *= method.direction * sign(det( [A; t'] ))
		t ./ norm(t)
	end

	function Delta(v, DHv = DH(v))
		[DHv; tang(DHv)'] \ [H(v); 0]
	end

	DHv = DH(method.v)
	v₀ = method.v + method.h * tang(DHv)
	DHv₀ = DH(v₀)

	dv = Delta(v₀, DHv₀)
	κ′ = norm(Delta(v₀ - dv)) / norm(dv)
	δ′ = norm(dv)
	α′ = acos( tang(DHv)' * tang(DH(v₀)) )[1]

	while norm(H(v₀)) > method.ɛ
		print(".") #
		v₁ = v₀ - Delta(v₀)
		v₀ = v₁
	end

	f = max(sqrt(κ′ / method.κ), sqrt(δ′ / method.δ), α′ / method.α)
	f = max(min(f, 2), 1/2)
	method.h = method.h / f;

	if (f < 2)
		method.v = v₀
	else
		print(" ") #
		m = max(sqrt(κ′ / method.κ), sqrt(δ′ / method.δ), α′ / method.α) #
		print(m == sqrt(κ′ / method.κ) ? "κ" : m == sqrt(δ′ / method.δ) ? "δ" : "α") #
		step(method, H, DH)
	end
end



# broyden

function initialize(method :: Broyden, H, DH, v₀)
	# TODO
end

function step(method :: Broyden, H, DH)
	# TODO
end



