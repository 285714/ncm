module NCM

abstract ContinuationMethod

type PredictorCorrector <: ContinuationMethod	# encapsulates the state of the Predictor-Corrector method
	h :: Real		# Stepsize
	eps :: Real		# Newton ϵ
	kNom :: Real	# κ
	dNom :: Real	# δ
	aNom :: Real	# α
	v :: Array		# current Position
	path :: Array   # v's history, should be saved in the generic cm method
end

PredictorCorrector(h, eps) = PredictorCorrector(h, eps, 0, 0, 0, [], [])


function initialize end
function intact end
function step end
function finish end

function cm(method :: ContinuationMethod,
            x,
			F :: Function,
			DF :: Function = jacobian(F),
			trace = x -> x)
	initialize(method, F, DF, x)

	n = 0;
	while intact(method, n)
		step(method, F, DF)
		n += 1
	end

	return finish(method, F, DF, x, n)
end


export PredictorCorrector, cm


function initialize(method :: PredictorCorrector, F, DF, x)
	method.v = x
	method.path = x
end

function intact(method :: PredictorCorrector, n)
	n < 1000
end

function step(method :: PredictorCorrector, F, DF)
	function tang(A)
		A2 = [A; rand(1,2)]
		t = A2 \ [zeros(1,1); 1]
		return  t ./ norm(t)
	end

	function Delta(v)
		[DF(v); tang(DF(v))'] \ [F(v); 0]
	end

	v0 = method.v + method.h * tang(DF(method.v))

	while norm(F(v0)) > method.eps
		v1 = v0 - Delta(v0)
		v0 = v1
	end

	method.v = v0
	method.path = [method.path v0];
end

function finish(method :: PredictorCorrector, F, DF, x, n)
	println("Finished path following")
	return method.path
end


end
