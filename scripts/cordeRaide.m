function [y] = cordeRaide(f,duration,L,Ft,rho,phi,kappa,sigma0,sigma1,xe,xrecord,fs)
% corde raide 
    k = 1/fs;
    Nt = floor(duration*fs);
    t = (0:Nt-1)*k;
    
    S = pi * phi.^2/4;           %surface de la section droite (m^2)
    
    gamma = sqrt(Ft./rho./S);
    assert(sum(isinf(gamma))==0,"gamma contient des valeurs inf");
    h = sqrt((max(gamma)^2 * k^2 + sqrt(max(gamma)^4 *k^4 + 16 * kappa^2 *k^2))/2)*1.05;
    

    Nx = floor(L/h);
    h = L/Nx;
    x = (0:Nx)*h;

    u_fut = zeros(size(x));
    u_pre = zeros(size(x));   
    u_pas = zeros(size(x));
    
    lambda = k*gamma/h;

    y = zeros(size(t));

    %boucle temporelle
    tic
    for n = 2:Nt-1
        %%%operateur d'interpolation
        l_e = floor(xe(n)/h)+1;           %partie entière
        alpha_e = xe(n)/h - (l_e-1);      %partie fractionnaire
        %spreading operator ordre 1
        J_l1 = zeros(size(x));
        J_l1(l_e) = (1-alpha_e)/h;
        J_l1(l_e+1) = alpha_e/h;
        %CL supportées
        l=2;
        u_fut(l) = (u_pre(l) * (2 - 2*lambda(n).^2 - 6* kappa^2*k^2/h^4 - 4*sigma1*k/h^2)...
            + u_pas(l)*(-1 + sigma0*k +4* sigma1 *k/h^2)...
            + u_pre(l+1)*(lambda(n).^2 + 4*kappa^2*k^2/h^4 + 2*sigma1*k/h^2)...
            -kappa^2*k^2/h^4*(u_pre(l+2)-u_pre(l))-2*sigma1*k/h^2*(u_pas(l+1))...
            + J_l1(l)*f(n)*k^2) / (1 + sigma0 * k);
        l=Nx;
        u_fut(l) = (u_pre(l) * (2 - 2*lambda(n).^2 - 6* kappa^2*k^2/h^4 - 4*sigma1*k/h^2)...
            + u_pas(l)*(-1 + sigma0*k +4* sigma1 *k/h^2)...
            +  u_pre(l-1) *(lambda(n).^2 + 4*kappa^2*k^2/h^4 + 2*sigma1*k/h^2)...
            -kappa^2*k^2/h^4*(-u_pre(l)+u_pre(l-2))-2*sigma1*k/h^2*u_pas(l-1)...
            + J_l1(l)*f(n)*k^2) / (1 + sigma0 * k);

        for l = 3:Nx-1
            u_fut(l) = (u_pre(l) * (2 - 2*lambda(n).^2 - 6* kappa^2*k^2/h^4 - 4*sigma1*k/h^2)...
                + u_pas(l)*(-1 + sigma0*k +4* sigma1 *k/h^2)...
                + (u_pre(l+1) + u_pre(l-1))*(lambda(n)^2 + 4*kappa^2*k^2/h^4 + 2*sigma1*k/h^2)...
                -kappa^2*k^2/h^4*(u_pre(l+2)+u_pre(l-2))-2*sigma1*k/h^2*(u_pas(l+1)+u_pas(l-1))...
                + J_l1(l)*f(n)*k^2) / (1 + sigma0 * k);
        end

        y(n) = u_pre(floor(xrecord(n)/h));
        u_pas = u_pre;
        u_pre = u_fut;
    end
end



