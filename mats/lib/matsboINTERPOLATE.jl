module matsboINTERPOLATE
export interpolateLanczos, interpolateTrigonometric

using matsboUTIL

# simple periodic (!) Lanczos interpolation
function interpolateLanczos(V,a::Integer)
	return matsboUTIL.vectorize() do y
		local sum = zero(V[1])
		local N = length(V)
		for i in floor(Integer,y)+(-a+1:a)
			local tmp = y-i
			sum += sinc(tmp)*sinc(tmp/a)*V[mod(i-1,N) + 1]
		end
		return sum
	end
end
interpolateLanczos(V, a::Integer, x) = interpolateLanczos(V,a)(x)



# returns trigonometric polynomial.
# use with 2a,-2b and divide by 2m+1 to use with rfft coefficients.
function interpolateTrigonometric(a₀, a, b)
	return matsboUTIL.vectorize() do x
		a₀ + sum(a.*cos(x*(1:length(a))) + b.*sin(x*(1:length(b))))
	end
end
interpolateTrigonometric(a₀, a, b, x) = interpolateTrigonometric(a₀, a, b)(x)

end #module
