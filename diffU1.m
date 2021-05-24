function Uval = diffU1(alpha,D1,N,i,j,ell,H,s,y)

% Derivative of U_1^(n)(x,y,s) [Equation (61)]
Uval = zeros(size(y));
for m = 0:N-1
    sig = mod(i,2); % Determine whether to include sinh or cosh
    lambda = m*pi/H;
    mu = sqrt(s/D1 + lambda^2); 
    if sig == 0
        Uval = Uval + (lambda)^j * alpha(m+1) * mu^i * sinh(mu*ell) * cos(lambda*y + j*pi/2);
    elseif sig == 1
        Uval = Uval + (lambda)^j * alpha(m+1) * mu^i * cosh(mu*ell) * cos(lambda*y + j*pi/2);
    end
end

end