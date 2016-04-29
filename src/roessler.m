function roessler(a, b, c, y0, T)

y = @(t,x) [ -x(2) - x(3);
              x(1) + a * x(2);
              b + x(3) * (x(1) - c) ];

[T,Y] = ode45(y, [0 T], y0);
hold on;
plot3(Y(:,1), Y(:,2), Y(:,3));
% scatter3(Y(end,1), Y(end,2), Y(end,3));

end