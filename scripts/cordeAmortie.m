function [y] = cordeAmortie(f,duration,L,Ft,rho,phi,sigma0,xe,xrecord,fs)
% corde amortie 
    k = 1/fs;
    Nt = floor(duration*fs);
    t = (0:Nt-1)*k;
    
    S = pi * phi.^2/4;           %surface de la section droite (m^2)
    gamma = sqrt(Ft./rho./S);
    
    %assert(sum(isinf(gamma))==0,"gamma contient des valeurs inf");
    h = max(gamma) * k; % condition de stabilité : k*gamma <= h

    Nx = floor(L/h);
    h = L/Nx;
    x = (0:Nx-1)*h;

    u_future = zeros(size(x));
    u_present = zeros(size(x));
    u_past = zeros(size(x));

    lambda = k*gamma/h;

    y = zeros(size(t));
    for n = 2:Nt-1
        % Spreading operator ordre 1
        J_ex = zeros(size(x));
        J_ex(floor(xe(n)/h)+1) = 1/h;

        % boucle temporelle
        for l = 2:Nx-1 
            u_future(l) = (u_present(l) * (2-2*lambda(n).^2) - u_past(l)*(1-sigma0*k) ...
                + lambda(n).^2*(u_present(l+1) + u_present(l-1) + J_ex(l)*f(n)*k^2)) / (1 + sigma0*k);
            
        end
        u_past = u_present;
        u_present = u_future;
        y(n) = u_future(floor(xrecord(n)/h)+1);
    end
end