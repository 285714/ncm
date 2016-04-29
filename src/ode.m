function ys = ode(f, y, h, N)

RK = [1/2 0 0 0; 0 1/2 0 0; 0 0 1 0; 1/6 2/6 2/6 1/6]; % classic
[s, ~] = size(RK);
[dim, ~] = size(y);
ys = NaN(dim,N);


for i=1:N
    ys(:,i) = y;
    k = zeros(dim, s);
    
    k(:,1) = f( y );
    for j=1:s-1
        k(:,j+1) = f( y + h * k * RK(j,:)' );
    end
    
    V = k * RK(s,:)';
    y = y + h*V;
end
