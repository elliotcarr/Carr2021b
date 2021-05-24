function U = higher_order_terms(N,w,wd,U1,U2,D1,D2,H,L,ell,x1,x2,y1,y2,M,s)

% Derivative of U_1^(0)(x,s) and U_2^(0)(x,s) [Equations (59) and (60)]
U01x = @(i) U1(ell,s,i);
U02x = @(i) U2(ell,s,i);

% Derivative of U_1^(n)(x,y,s) and U_2^(n)(x,y,s) [Equations (61) and (62)]
diffU1x = @(alpha,i,j,y) diffU1(alpha,D1,M,i,j,ell,H,s,y);
diffU2x = @(beta,i,j,y) diffU2(beta,D2,M,i,j,ell,H,L,s,y);

% F^(n)(x,s) and G^(n)(x,s) [Equations (36) and (37)]
F = @(alpha,beta,m,y) Ffunc(alpha,beta,w,diffU1x,diffU2x,U01x,U02x,m,y);
G = @(alpha,beta,m,y) Gfunc(alpha,beta,D1,D2,w,wd,diffU1x,diffU2x,U01x,U02x,m,y);

U = zeros(size([x1;x2],1),N-1);
alpha = zeros(N-1,M);
beta = zeros(N-1,M);
shift = length(x1);
for n = 1:N-1

    for m = 0:M-1
        
        % lambda_m, mu_{1,m}(s) and mu_{2,m}(s) defined on page 11
        lambda = m*pi/H;
        mu1 = sqrt(s/D1 + lambda^2);
        mu2 = sqrt(s/D2 + lambda^2);
        
        % \tilde{F}_{m}^(n)(s) and \tilde{G}_{m}^(n)(s) defined on page 8
        Ft = integral(@(y) F(alpha,beta,n,y) .* cos(lambda*y),0,H); 
        Gt = integral(@(y) G(alpha,beta,n,y) .* cos(lambda*y),0,H);
        
        % \tilde{alpha}_m^(n)(s) and \tilde{beta}_m^(n)(s) defined on page 11
        if m > 0
            alphat = 2/(H*D1*mu1*cosh(mu1*ell));
            betat = 2/(H*D2*mu2*sinh(mu2*(ell-L)));
        
            numratr = (2*Ft/H) - betat*Gt*cosh(mu2*(ell-L));
            denmntr = alphat*sinh(mu1*ell) - betat*cosh(mu2*(ell-L));
        else
            alphat = 1/(H*D1*mu1*cosh(mu1*ell));
            betat = 1/(H*D2*mu2*sinh(mu2*(ell-L)));
        
            numratr = (Ft/H) - betat*Gt*cosh(mu2*(ell-L));
            denmntr = alphat*sinh(mu1*ell) - betat*cosh(mu2*(ell-L));
        end
        
        V = numratr/denmntr; % V_m^(n)(s) defined on page 8
        
        % alpha_m^(n)(s) and beta_m^(n)(s) defined on page 8
        alpha(n,m+1) = alphat*V; 
        beta(n,m+1) = betat*(V - Gt);
        
        % U_1^(n)(x,y,s) and U_2^(n)(x,y,s) [Equations (55) and (56)]
        U(1:shift,n) = U(1:shift,n) + alpha(n,m+1)*sinh(mu1*x1).*cos(lambda*y1);
        U(shift+1:end,n) = U(shift+1:end,n) + beta(n,m+1)*cosh(mu2*(x2-L)).*cos(lambda*y2);
    end
end
end
