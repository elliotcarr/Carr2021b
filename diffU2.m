function Uval = diffU2(beta,D2,N,i,j,ell,H,L,s,y)

% Derivative of U_2^(n)(x,s) [Equation (62)]
Uval = zeros(size(y)); 
for m = 0:N-1
    sig = mod(i,2); % Determine wheher to include sinh or cosh
    lambda = m*pi/H;
    mu = sqrt(s/D2 + lambda^2);    
    if sig == 0
        Uval = Uval + (lambda)^j * beta(m+1) * mu^i * cosh(mu*(ell-L)) * cos(lambda*y + j*pi/2); 
    elseif sig == 1
        Uval = Uval + (lambda)^j * beta(m+1) * mu^i * sinh(mu*(ell-L)) * cos(lambda*y + j*pi/2); 
    end  
end

end