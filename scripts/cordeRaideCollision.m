function [out] = cordeRaideCollision(f,duration,L,Ft,rho,diametre,kappa,sigma0,sigma1,xe,xc,xrecord,fs)
% corde raide avec collision à un ressort au point xc
    k = 1/fs;
    Nt = floor(duration*fs);
    t = (0:Nt-1)*k;
    
    S = pi * diametre.^2/4;           %surface de la section droite (m^2)
    
    gamma = sqrt(Ft./rho./S);
    assert(sum(isinf(gamma))==0,"gamma contient des valeurs inf");
    h = sqrt((max(gamma)^2 * k^2 + sqrt(max(gamma)^4 *k^4 + 16 * kappa^2 *k^2))/2)*1.05;
    
    alpha = 1; % ordre du ressort collision

    Nx = floor(L/h);
    h = L/Nx;
    x = (0:Nx)*h;

    u_fut = zeros(size(x));
    u_pre = [(1:floor(Nx/2)) (Nx:-1:floor(Nx/2))-floor(Nx/2)]*(-0.5);
    u_pas = [(1:floor(Nx/2)) (Nx:-1:floor(Nx/2))-floor(Nx/2)]*(-0.5);
    

    
    lambda = k*gamma/h;
    eps = 10^(-6);
    yc = 1 ; % position d'équilibre du ressort collision
    kc = 8;
    fc = zeros(size(t));
    out = zeros(size(t));
    
    eta_fut = yc;
    %boucle temporelle
    tic
    for n = 2:Nt-1
        %%%operateur d'interpolation
        l_e = floor(xe(n)/h)+1;
        l_c = floor(xc(n)/h)+1;
        alpha_e = xe(n)/h - (l_e-1);
        alpha_c = xc(n)/h - (l_c-1);
        %spreading operator ordre 0
        J_l1 = zeros(size(x));
        J_l1(l_e) = (1-alpha_e)/h;
        J_l1(l_e+1) = alpha_e/h;

        J_l0 = zeros(size(x));
        J_l0(l_c) = 1/h;
        %
        
        eta_pas = (1- alpha_c) * u_pas(l_c) + alpha_c * u_pas(l_c+1);
        eta_pre = (1- alpha_c) * u_pre(l_c) + alpha_c * u_pre(l_c+1);
        eta_pas_lp1 = (1- alpha_c) * u_pas(l_c+1) + alpha_c * u_pas(l_c+2);
        eta_pre_lp1 = (1- alpha_c) * u_pre(l_c+1) + alpha_c * u_pre(l_c+2);
        eta_pas_lm1 = (1- alpha_c) * u_pas(l_c-1) + alpha_c * u_pas(l_c);
        eta_pre_lm1 = (1- alpha_c) * u_pre(l_c-1) + alpha_c * u_pre(l_c);
        eta_pre_lp2 = (1- alpha_c) * u_pre(l_c+2) + alpha_c * u_pre(l_c+3);
        eta_pre_lm2 = (1- alpha_c) * u_pre(l_c-2) + alpha_c * u_pre(l_c-1);
        
        b = (eta_pre*(-2 + 2*lambda(n).^2+6*kappa^2*k^2/h^4+4*sigma1*k/h^2) ...
            + eta_pas * (2 - sigma0*k-4*sigma1*k/h^2) ...
            - (eta_pre_lp1 + eta_pre_lm1) * (lambda(n)^2 + 4*kappa^2*k^2/h^4+2*sigma1*k/h^2) ...
            + (eta_pas_lp1 + eta_pas_lm1) * 2*sigma1*k/h^2 ...
            + (eta_pre_lp2 + eta_pre_lm2) * k^2*kappa^2/h^4);
        

        
        % Newton-Raphson
        % phi_a = kc * positivePart(eta_pas - yc)^(alpha+1)/(alpha+1);
        m = k^2/rho(n)/S(n)/(1+sigma0*k);
        n_iter = 0;
        dif = 1;
        a = eta_pas;
        r = eta_fut - a; % estimation Raphson
        while dif>eps
            %phi_ra = kc * positivePart(eta_fut - yc)^(alpha+1)/(alpha+1);
            if abs(r) < 10^(-12) % approximation pour eviter les instabilités numériques quand on divise par r
                G = r + m*kc*(a - yc)^alpha + b;
                dG = 1 - m*alpha*kc*(a-yc)^(alpha-1);
            else
                G = r + m / r * (phi(r+a,yc,kc,alpha) - phi(a,yc,kc,alpha)) + b;
                dG = 1 - m/r^2*(phi(r+a,yc,kc,alpha) - phi(a,yc,kc,alpha)) + m*kc/r*(r + a - yc)^alpha;
            end
            r_new = r - G/dG;
            dif = abs(r_new-r);
            r = r_new;
            n_iter = n_iter + 1;
            if n_iter>200
                error('nombre d iteration trop important')
            end
        end
        eta_fut = r + a;
        fc(n) = (positivePart(eta_fut - yc)^(alpha+1) - positivePart(eta_pas - yc)^(alpha+1)) * kc /(alpha+1)/r;

        %CL supportées
        l=2;
        u_fut(l) = (u_pre(l) * (2 - 2*lambda(n).^2 - 6* kappa^2*k^2/h^4 - 4*sigma1*k/h^2)...
            + u_pas(l)*(-1 + sigma0*k +4* sigma1 *k/h^2)...
            + u_pre(l+1)*(lambda(n).^2 + 4*kappa^2*k^2/h^4 + 2*sigma1*k/h^2)...
            -kappa^2*k^2/h^4*(u_pre(l+2)-u_pre(l))-2*sigma1*k/h^2*(u_pas(l+1))...
            + J_l1(l)*f(n)*k^2  + J_l0(l)*fc(n)*k^2/rho(n)/S(n)) / (1 + sigma0 * k);
        l=Nx;
        u_fut(l) = (u_pre(l) * (2 - 2*lambda(n).^2 - 6* kappa^2*k^2/h^4 - 4*sigma1*k/h^2)...
            + u_pas(l)*(-1 + sigma0*k +4* sigma1 *k/h^2)...
            +  u_pre(l-1) *(lambda(n).^2 + 4*kappa^2*k^2/h^4 + 2*sigma1*k/h^2)...
            -kappa^2*k^2/h^4*(-u_pre(l)+u_pre(l-2))-2*sigma1*k/h^2*u_pas(l-1)...
            + J_l1(l)*f(n)*k^2  + J_l0(l)*fc(n)*k^2/rho(n)/S(n)) / (1 + sigma0 * k);

        for l = 3:Nx-1
            u_fut(l) = (u_pre(l) * (2 - 2*lambda(n).^2 - 6* kappa^2*k^2/h^4 - 4*sigma1*k/h^2)...
                + u_pas(l)*(-1 + sigma0*k +4* sigma1 *k/h^2)...
                + (u_pre(l+1) + u_pre(l-1))*(lambda(n)^2 + 4*kappa^2*k^2/h^4 + 2*sigma1*k/h^2)...
                -kappa^2*k^2/h^4*(u_pre(l+2)+u_pre(l-2))-2*sigma1*k/h^2*(u_pas(l+1)+u_pas(l-1))...
                + J_l1(l)*f(n)*k^2 + J_l0(l)*fc(n)*k^2/rho(n)/S(n)) / (1 + sigma0 * k);
        end

        out(n) = u_pre(floor(xrecord(n)/h));
        u_pas = u_pre;
        u_pre = u_fut;
    end
    yyaxis left
    plot(t,out);
    yyaxis right 
    plot(t,fc);
end

function [w] = phi(u, yc, K, alpha)
    v = u - yc;
    if v <= 0 % partie positive
        v = 0;
    end
    w = K*v^(alpha+1)/(alpha+1);
end