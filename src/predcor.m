function [us,u0] = predcor(F, DF, accept, u0)

h = 0.01;
newtonEps = 0.01;
N = length(u0) - 1;
us = [];

kNom = 1.5; %1.0
dNom = 1.5; %0.5
aNom = 1.5; %0.6


warning('off', 'MATLAB:warn_truncate_for_loop_index');
for i = 0:Inf
    us = [us u0];
    if (accept(i, u0))
        break;
    end
    
    % Predictor
    
    v0 = u0 + h * t(DF(u0));        % Euler
    % plot([u0(end), v0(end)], [norm(u0(1:N)), norm(v0(1:N))], '--', 'Color', [0,0,0]);
    
    dv = Delta(v0);
    kappa = norm(Delta(v0 - h * dv)) / norm(dv);
    delta = norm(dv);
    alpha = acos( t(DF(u0))' * t(DF(v0)) );
    
    
%     % Corrector
    
    while norm(F(v0)) > newtonEps
        dx = Delta(v0);
        v1 = v0 - dx;
        v0 = v1;
        
        fprintf('>');
    end
    
    
    % adapt step size
    
    f = max([ sqrt(kappa/kNom), sqrt(delta/dNom), alpha/aNom ]);
    f = max(min(f,2),1/2);
    h = h / f;
    
    if (f < 2)
        u0 = v0;
    end
end



    function dv = Delta(v)
        dv = linsolve([DF(v); t(DF(v))'], [F(v); 0]);
    end

    function t = t(A)
        if rank(A) < N    % singular point..
            fprintf('#');
            t = rand(N+1, 1) * newtonEps;
            return;
        end

        while true
            A2 = [A; rand(1,N+1)];

            if (rank(A2) == N+1)
                break;
            end
        end

        t = linsolve(A2, [zeros(N,1); 1]);
        t = t ./ norm(t);

        if (det( [A; t'] ) < 0)
            t = -t;
        end
    end

end