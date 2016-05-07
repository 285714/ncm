
#	convert input/output vector
# function rToC(V)
#   local N = (length(V)-1)÷3
#   local c(x) = [ x[1]+0im; complex(x[2:end÷2+1], x[end÷2+2:end]) ]
#   return matsboUTIL.inputConvert(V, [1;N+1;2N+1;3N+1;3N+2], [c,c,c,identity])
# end
#
# function cToR(V...)
#   local d(x) = [real(x[1]); real(x[2:end]); imag(x[2:end])]
#   return matsboUTIL.outputConvert([d,d,d,identity], V...)
# end

# optimization target for roessler
#		operating on real Fourier coefficients plus ω
# function Hroessler(V::Vector{Float64})
#   local X,Y,Z,ω
#   X,Y,Z,ω = rToC(V)
#   local N = (length(X)-1)*2
#
#   local A,B,C
# 	local D = im*ω.*(0:N÷2) #deriv
#   A = -Y - Z - D.*X
#   B = X + a*Y - D.*Y
#   C = [N*b; zeros(N÷2)] - c*Z + rfft(irfft(Z,N).*irfft(X,N)) - D.*Z
#
#   return cToR(A,B,C,reduce((X,x)->X+real(x), .0, X))
# end

# derivative
#Jroessler(v) = matsboNWTN.forwardDifference(Hroessler, v, ϵ=1e-4)


# generic objective function idea...
# function Htest(V)
# 	local X,Y,Z,ω
# 	X,Y,Z,ω = rToC(V)
# 	local N = (length(X)-1)*2
#
# 	local tmp = mapslices(roessler, irfft([X Y Z], N, [1]), [2])
# 	local tmp2 = tmp[1,1] # not working?! why???
# 	tmp = rfft(tmp, [1])
#
# 	local A,B,C
# 	local D = im*ω.*(0:N÷2)
# 	A = tmp[:,1]-D.*X
# 	B = tmp[:,2]-D.*Y
# 	C = tmp[:,3]-D.*Z
#
# 	return cToR(A,B,C,tmp2)
# end
#
# Jtest(v) = matsboNWTN.forwardDifference(Htest, v, ϵ=1e-4)
