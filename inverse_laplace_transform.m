function invL = inverse_laplace_transform(G,t,Np,c,z)

% Computes inverse Laplace transform L^(-1){G(s)} for some specified function
% G(s) based on Trefethen et al. (2006)
invL = 0.0;
for k = 1:2:Np
    s = z(k)/t;
    invL = invL-c(k)*G(s);
end
invL = 2*real(invL)/t;

end
