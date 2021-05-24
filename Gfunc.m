function Gval = Gfunc(alpha,beta,D1,D2,w,gp,Uj1x,Uj2x,U01x,U02x,m,y)

% G^(n)(x,s) [Equation (37)]
Gval = (w(y).^m/factorial(m)).*(D2*U02x(m+1) - D1*U01x(m+1));
for k = 1:m-1
    Gval = Gval + (w(y).^k./factorial(k)).*(D2*Uj2x(beta(m-k,:),k+1,0,y) ...
        - D1*Uj1x(alpha(m-k,:),k+1,0,y)) ...
        - gp(y).*(w(y).^(k-1)/factorial(k-1)).*(D2*Uj2x(beta(m-k,:),k-1,1,y) ...
        - D1*Uj1x(alpha(m-k,:),k-1,1,y));
end

end