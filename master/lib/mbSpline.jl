type Spline
  d
  U #knots
  C #coefficients
  function Spline(d,U,C)
    p = sortperm(U)
    U,C = U[p], C[p]
    U = [ fill(U[1],d); U; fill(U[end],d) ]
    C = [ fill(C[1],d); C; fill(C[end],d) ]
    return new(d,U,C)
  end
end

test = Spline(1,linspace(0,1,10),rand(10))

function fastfind(A,x)
  i,j = 1,length(A)
  A[i] > x && return 0
  A[j] < x && return j
  while abs(i-j) > 1
    k = (i+j) ÷ 2
    if A[k] < x; i = k
    else j = k end
  end
  return i
end

function splinematrix(S::Spline,i,k,x)
  d,U = S.d,S.U
  M = zeros(k,k+1)
  for a in 1:k
    M[a,a] = (U[i+a] - x) / (U[i+a] - U[i+a-k])
    M[a,a+1] = (x - U[i+a-k]) / (U[i+a] - U[i+a-k])
  end
  return M
end

splinematrix(test,2,1,.1)

@noinline function call(S::Spline, x)
  core(S::Spline, x)
end

@noinline function core(S::Spline, x)
  d,U,C = S.d,S.U,S.C
  i = fastfind(U,x)
  # i = clamp(i,d+1,length(U)-d-1)
  i = max(i,d+1)
  println("x")
  c = C[i-d:i]
  for k in d:-1:1
    c = splinematrix(S,i,k,x) * c
  end
  return c[1]
end

using PyPlot
plot(linspace(0,1,length(test.C)), test.C)
plot(linspace(0,1,100), map(x->core(test, x), linspace(test.U[1],test.U[end],100)))

splinematrix(test, 10, 1, 1.0)

x = sort(rand(10))












fastfind(x,.0)
