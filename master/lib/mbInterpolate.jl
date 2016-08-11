__precompile__()

module mbInterpolate
export interpolateLanczos, interpolateTrigonometric

using mbUtil

"""
    interpolateLanczos(V, a::Integer)

simple periodic (!) Lanczos interpolation
"""
function interpolateLanczos(V,a::Integer)
	return mbUtil.vectorize() do y
		sum = zero(V[1])
		N = length(V)
		for i in floor(Integer,y)+(-a+1:a)
			local tmp = y-i
			sum += sinc(tmp)*sinc(tmp/a)*V[mod(i-1,N) + 1]
		end
		return sum
	end
end
interpolateLanczos(V, a::Integer, x) = interpolateLanczos(V,a)(x)



"""
    interpolateTrigonometric(a₀, a, b)
Returns trigonometric polynomial.
Use with 2a,-2b and divide by 2m+1 to use with rfft coefficients.
"""
function interpolateTrigonometric(a₀, a, b)
	return mbUtil.vectorize() do x
		a₀ + reduce(0, 1:length(a)) do I,i
			tmp = i*x
			I + a[i]*cos(tmp) + b[i]*sin(tmp)
		end
	end
end
interpolateTrigonometric(a₀, a, b, x) = interpolateTrigonometric(a₀, a, b)(x)

end #module
