module matsboRK
export rk1, rk4

import Base.length

type ButcherTableau
	A::Array{Float64,2}
	b::Array{Float64,1}
	c::Array{Float64,1}
end

length(BT::ButcherTableau) = length(BT.b)

function rk(BT::ButcherTableau)
	return function rk(f::Function, t₀::Float64, y₀::Vector{Float64}, h::Float64, pred::Function;
		init=Void, callback=Void)

		local A,b,c,s,t,y,K
		A,b,c = BT.A, BT.b, BT.c
		s = length(BT)
		K = cell(s)

		t,y = deepcopy(t₀), deepcopy(y₀)
		K[1] = f(t,y)

		init≠Void && init()

		while pred(t,y,K[1])
			for i in 2:s
				K[i] = f( t + h*c[i], y + h*(A[i,1:i-1]*K[1:i-1])[1] )
			end

			t += h
			y += h*sum(map(*, K, b))
			K[1] = f(t,y)

			callback≠Void && callback(t,y,K)
		end

		return Void
	end
end




BTrk1 = ButcherTableau(
	zeros(1,1),
	ones(1),
	zeros(1)
)
rk1 = rk(BTrk1)

BTrk4 = ButcherTableau(
	[
		.0	.0	.0	.0
		.5	.0	.0	.0
		.0	.5	.0	.0
		.0	.0	1.0	.0
	],
	[1/6; 1/3; 1/3; 1/6],
	[.0; .5; .5; 1.0]
)
rk4 = rk(BTrk4)


end
