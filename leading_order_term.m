function [U1,U2] = leading_order_term(D1,D2,ell,L,C0,QL)

% mu_1(s) and mu_2(s) defined on page 11
mu1 = @(s) sqrt(s/D1);
mu2 = @(s) sqrt(s/D2);

% psi_1(s) and psi_2(s) defined on page 11
psi1 = @(s) -D1*mu1(s)*(exp(-mu1(s)*ell) + exp(mu1(s)*ell));
psi2 = @(s) -D2*mu2(s)*(exp(-mu2(s)*(L-ell)) - exp(mu2(s)*(L - ell)));

% gamma_11(s),...,gamma_24(s) defined on page 11 
gamma11 = @(s) -D1*mu1(s)*exp(-mu1(s)*ell)*C0(s)/psi1(s);
gamma12 = @(s) -1/psi1(s);
gamma13 = @(s) -D1*mu1(s)*exp(mu1(s)*ell)*C0(s)/psi1(s);
gamma14 = @(s) 1/psi1(s);
gamma21 = @(s) D2*mu2(s)*exp(-mu2(s)*ell)*QL(s)/psi2(s);
gamma22 = @(s) -exp(-mu2(s)*L)/psi2(s);
gamma24 = @(s) -exp(mu2(s)*L)/psi2(s);
gamma23 = @(s) D2*mu2(s)*exp(mu2(s)*ell)*QL(s)/psi2(s);

% V^(0)(s) defined on page 7
V0 = @(s) (gamma11(s)*exp(mu1(s)*ell) + gamma13(s)*exp(-mu1(s)*ell) ...
    - gamma21(s)*exp(mu2(s)*ell) - gamma23(s)*exp(-mu2(s)*ell)) / ...
    (gamma22(s)*exp(mu2(s)*ell) + gamma24(s)*exp(-mu2(s)*ell) ...
    - gamma12(s)*exp(mu1(s)*ell) - gamma14(s)*exp(-mu1(s)*ell));

% U_1^(0)(x,s) and U_2^(0)(x,s) [Equations (47) and (48)] (i=0) and 
% derivative of U_1^(0)(x,s) and U_2^(0)(x,s) [Equations (59) and (60)]
U1 = @(x,s,i) (gamma11(s) + gamma12(s)*V0(s))*(mu1(s)^i)*exp(mu1(s)*x) + (gamma13(s) + gamma14(s)*V0(s))*(-mu1(s))^i*exp(-mu1(s)*x);
U2 = @(x,s,i) (gamma21(s) + gamma22(s)*V0(s))*(mu2(s)^i)*exp(mu2(s)*x) + (gamma23(s) + gamma24(s)*V0(s))*(-mu2(s))^i*exp(-mu2(s)*x);

end