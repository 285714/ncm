
using PyPlot
using ODE
using Debug


const m = 3;

const a = 0.1;
const b = 0.1;
const c = 4;


function splitC(C)
	ω   = C[ 1 ]
	cxr = C[ 2     :  2*m+1]
	cxr = [cxr; -sum(cxr)]
	cxi = C[ 2*m+2 :  4*m+2]
	cyr = C[ 4*m+3 :  6*m+3]
	cyi = C[ 6*m+4 :  8*m+4]
	czr = C[ 8*m+5 : 10*m+5]
	czi = C[10*m+6 : 12*m+6]
	
	return (ω,
	        cxr + im*cxi,
	        cyr + im*cyi,
	        czr + im*czi)
end

function f(C)
	ω, cx, cy, cz = splitC(C)

	k         = -m:m
	sumMatrix = reverse(cx) * transpose(cz)
	s_         = map(d -> sum(diag(sumMatrix, d)), k)

	eq1 = ω .* cx .* im .* k  +  cy + cz						# (1)
	eq2 = ω .* cy .* im .* k  -  cx - a .* cy;					# (2)
	eq3 = ω .* cz .* im .* k  +  b .* (k.==0) - c .* cz  -  s_	# (3)

	return [real(eq1); imag(eq1); real(eq2); imag(eq2); real(eq3); imag(eq3)]
end

function Df(C)
	l = length(C)
	D = zeros(l, l);
	eps = 0.1;

	for i = 1:l
		D[:,i] = f(C + eps .* ((1:l).==i)) ./ eps;
	end

	return D
end

function newton(C0, eps)
	while norm(f(C0)) > eps
		dC = Df(C0) \ f(C0);
		C0 = C0 - dC
	end

	return C0
end

function four(c, s)
	k = -m:m
	return exp(im .* s * k') * c
end



# Define time vector and interval grid
const dt = 0.001
const tf = 6.0
t = 0:dt:tf

# Initial position in space
# r0 = [0.; 0.2; 0.]
# r0 = [-7.707; 0.451; 0.007336]
r0 = [-6.0; -0.243; 0]

# ode-solution:
R(t, x) = [-x[2] - x[3] ; x[1] + a*x[2] ; b + x[3]*(x[1]-c)]

t, pos = ode45(R, r0, t)
x = map(v -> v[1], pos)
y = map(v -> v[2], pos)
z = map(v -> v[3], pos)

plot3D(x, y, z)
show()


# initial guess:
s = 0:0.1:6.3
plot3D(four([0;0;0;0;1+10im;0;0],s), four([0;0;0;0;7+5im;0;0],s), four([0;0;0;0;0;0;0],s))


C0 = rand(12*m+6)
#=C0 = [3.1; 
      0;0;0;0;1;0;
      0;0;0;0;10;0;0;
      0;0;0;0;7;0;0;
      0;0;0;0;5;0;0;
      0;0;0;0;0;0;0;
      0;0;0;0;0;0;0]=#
print("-->")
C  = newton(C0, 0.00001)
println("<--")
ω, cx, cy, cz = splitC(C)

println(ω)
println(cx)
println(cy)
println(cz)

s = 0:ω/10:2*ω
plot3D(four(cx, s), four(cy, s), four(cz, s))

