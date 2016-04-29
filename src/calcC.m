function C = calcC(a, b, c)

m = 20;


% y0 = [-0.08438;
%       -6.566;
%        0.01944];
% T = 5.9
h = 0.1;
% T = 0:h:6;

y0 = [0.003275;
      -9.675;
      0.01396];
T = 0:h:12;

y = @(x) [ -x(2) - x(3);
            x(1) + a * x(2);
            b + x(3) * (x(1) - c) ];

% hold on; 
Y = ode(y, y0, h, length(T));
Y1 = Y(1,:);
Y2 = Y(2,:);
Y3 = Y(3,:);

plot3(Y1, Y2, Y3);
% scatter3(Y(end,1), Y(end,2), Y(end,3));



    function x = four(C, s)
        m = (length(C) - 1) / 2;
        k = (-m:m);
        x = exp(1i .* s * k) * C;
    end





    function c = fcoef(Y, T)
        c = NaN(2*m + 1, 1);

        for k = -m:m
            w = 2*pi / T(end);
            c(k+m+1,1) = 0.05 .* 1./(2*pi) .* Y * exp(-1i .* k * T' .* w);
        end
    end



cx = fcoef(Y1, T);
cy = fcoef(Y2, T);
cz = fcoef(Y3, T);

x = four(cx, T');
y = four(cy, T');
z = four(cz, T');

hold on;
 plot3(x, y, z);
 legend('ode', 'approx');

%  figure;
%  hold on;
%  plot(T, Y1, '-x');
% plot(T, x, '-x');
%  plot(T, Y2, '-x');
%  plot(T, y, '-x');
%  plot(T, Y3, '-x');
%  plot(T, z, '-x');
%  legend('1:ode', '1:approx', '2:ode', '2:approx', '3:ode', '3:approx');


cxr = real(cx);
cxi = imag(cx);
cyr = real(cy);
cyi = imag(cy);
czr = real(cz);
czi = imag(cz);

C = [12;
     cxr(1:end-1);
     cxi;
     cyr;
     cyi;
     czr;
     czi];

end