function [x,w,P]=gpopsGetLGR(N)

N = N-1;
N1 = N+1;
% Initial guess for LGR nodes
x = -cos(2*pi*(0:N)/(2*N+1))';

% The Legendre Vandermonde Matrix
P = zeros(N1,N1+1);
xold = 2;

% Free abscissae
free = 2:N1;

while max(abs(x-xold))>eps
    xold = x;
    P(1,:) = (-1).^(0:N1);
    P(free,1) = 1;
    P(free,2) = x(free);
    for k = 2:N1
        P(free,k+1) = ( (2*k-1)*x(free).*P(free,k)-(k-1)*P(free,k-1) )/k;
    end
    x(free) = xold(free)-((1-xold(free))/N1).*(P(free,N1)+P(free,N1+1))./(P(free,N1)-P(free,N1+1));
end

% The Legendre-Gauss-Radau Vandermonde
P = P(1:N1,1:N1);

% Compute the weights
w = zeros(N1,1);
w(1) = 2/N1^2;
w(free) = (1-x(free))./(N1*P(free,N1)).^2;

