function Fval = Ffunc(alpha,beta,w,Uj1x,Uj2x,U01x,U02x,m,y)

% F^(n)(x,s) [Equation (36)]
Fval = ((w(y).^m)./factorial(m)).*(U02x(m) - U01x(m));
for k = 1:m-1
    Fval = Fval + ((w(y).^k)./factorial(k)).*(Uj2x(beta(m-k,:),k,0,y) ...
        - Uj1x(alpha(m-k,:),k,0,y));
end

end