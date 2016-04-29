function galerkin

a = 0.1;
b = 0.1;

% m = 10;
% c = 4;
% T = 6;

m = 3;
c = 6.1;

    function f = f(C)
        omega = C( 1 );
        cxr = C( 2     :  2*m+1);
        cxr = [cxr; -sum(cxr)];
        cxi = C( 2*m+2 :  4*m+2);
        cyr = C( 4*m+3 :  6*m+3);
        cyi = C( 6*m+4 :  8*m+4);
        czr = C( 8*m+5 : 10*m+5);
        czi = C(10*m+6 : 12*m+6);
        cx  = cxr + 1i * cxi;
        cy  = cyr + 1i * cyi;
        cz  = czr + 1i * czi;

%         cx = C(1 : 2*m+1, :);
%         cy = C(2*m+2 : 4*m+2, :);
%         cz = C(4*m+3 : 6*m+3, :);

        k = (-m:m)';

%         s2 = zeros(2*m+1, 1);
%         for i = k'
%             for j = k'
%                 if abs(i+j) > m
%                     continue;
%                 end
%                 
%                 s2(m+1 + i+j, 1) = s2(m+1 + i+j, 1) + ...
%                     cx(m+1 + i,1) * cz(m+1 + j,1);
%             end
%         end

        sumMatrix = fliplr(cx')' * transpose(cz);
        s = NaN(2*m+1, 1);
        for d = k'
            s(m+1+d,1) = sum(diag(sumMatrix, d));
        end

        eq1 = omega .* cx .* 1i .* k  +  cy  +  cz;                      % (1)
        eq2 = omega .* cy .* 1i .* k  -  cx  -  a .* cy;                 % (2)
        eq3 = omega .* cz .* 1i .* k  -  b .* (k==0)  +  c .* cz  -  s;  % (3)
        
        f = [real(eq1); imag(eq1); real(eq2); imag(eq2); real(eq3); imag(eq3)];
    end

    function Df = Df(C)
        Df = DF(C, @f);
    end

    function DF = DF(C, f)
        fC = f(C);
        h = length(fC);
        w = length(C);
        DF = NaN(h, w);
        eps = 0.1;
        
        for i = 1:w
            DF(:,i) = ( f(C + eps .* ((1:w)'==i)) - fC ) ./ eps;
        end
    end
        

    function x0 = newton(x0, eps)
        while norm(f(x0)) > eps
            dx = Df(x0) \ f(x0);
            x0 = x0 - dx;
            
            disp(norm(f(x0)));
            showC(x0);
            hold off;
            drawnow;
        end
    end



% determine initial guess..
C0 = calcC(a, b, c);


% C0 = [3.1; 
%       0;0;0;0;1;0;
%       0;0;0;0;10;0;0;
%       0;0;0;0;7;0;0;
%       0;0;0;0;5;0;0;
%       0;0;0;0;0;0;0;
%       0;0;0;0;0;0;0];
% 
% C0 = 1000 * rand(12*m+6,1);
% C0 = C0 + 0.01 * rand(12*m+6,1);


    function x = four(C, s)
        k = (-m:m);
        x = exp(1i .* s * k) * C;
    end

    function showC(C)
        omega = C( 1 );
        cxr = C( 2     :  2*m+1);
        cxr = [cxr; -sum(cxr)];
        cxi = C( 2*m+2 :  4*m+2);
        cyr = C( 4*m+3 :  6*m+3);
        cyi = C( 6*m+4 :  8*m+4);
        czr = C( 8*m+5 : 10*m+5);
        czi = C(10*m+6 : 12*m+6);
        cx  = cxr + 1i * cxi;
        cy  = cyr + 1i * cyi;
        cz  = czr + 1i * czi;
        
        sgn = sign(omega);
        omega = abs(omega);
        s = sgn .* (0:0.01:2*omega)';
        x = four(cx, s);
        y = four(cy, s);
        z = four(cz, s);
        % scale = 1 ./ max([x; y; z]);
        
        warning('off', 'MATLAB:plot:IgnoreImaginaryXYZPart');
        plot3(x, y, z);
        hold on;
        % warning('off', 'MATLAB:specgraph:private:specgraph:UsingOnlyRealComponentOfComplexData');
        % scatter3([x(1) x(end)], [y(1) y(end)], [z(1) z(end)]);
        
        s_ = sgn .* (2*omega:0.1:20*omega)';
        x_ = four(cx, s_);
        y_ = four(cy, s_);
        z_ = four(cz, s_);
        warning('off', 'MATLAB:plot:IgnoreImaginaryXYZPart');
        plot3(x_, y_, z_, '--');
    end


for m_ = 20:20
    disp(['###### m = ', int2str(m), ' ##########']);
    m = m_
    tic;
    C = newton(C0, 2.0);
    % C = newton(C0, 1.8 / m^2);
    toc;
    disp(norm(f(C)));
    
    
    omega = C( 1 );
    cxr = C( 2     :  2*m+1);
    cxi = C( 2*m+2 :  4*m+2);
    cyr = C( 4*m+3 :  6*m+3);
    cyi = C( 6*m+4 :  8*m+4);
    czr = C( 8*m+5 : 10*m+5);
    czi = C(10*m+6 : 12*m+6);
    
    C1 = [omega;
          0; cxr; 0;
          0; cxi; 0;
          0; cyr; 0;
          0; cyi; 0;
          0; czr; 0;
          0; czi; 0];

%     figure;
%     hold on;
%     showC(C);
%     showC(C0);
%     
%     y0 = [6.813; 0.4996; 0.9644];
%     T = 2 * 3.1;
%     roessler(a, b, c, y0, T);
% 
%     legend('approx', '-', 'initial guess', '-', 'ode');
    
    C0 = C1;
end

figure;
hold on;
showC(C);

% y0 = [6.813; 0.4996; 0.9644];
y0 = [0.003275;
      -9.675;
      0.01396];
T = 4 * 3.0;
roessler(a, b, c, y0, T);

legend('approx', '-', 'ode');

disp([cxr+1i*cxi, cyr+1i*cyi, + czr+1i*czi]);

% compare with initial guess
% m = 3;
% showC(calcC(a,b,c));
% m=20;




    function g = g(u)
        c = u(1,1);
        g = f(u(2:end,1));
    end

    function Dg = Dg(u)
        Dg = DF(u, @g);
    end

    function a = accept(i,u)
        disp(i);
        disp(u(1,1));
        a = u(1,1) > 7;
        
        showC(u(2:end));
        title(u(1,1));
        drawnow;
        hold off;
    end

% [us,u0] = predcor(@g, ...
%                   @Dg, ...
%                   @accept, ...
%                   [c; C]);
% 
% for i = 1:size(us,2);
%     showC(us(2:end,i));
%     title(us(1,i));
%     hold off;
% end


end