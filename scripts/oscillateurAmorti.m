function [u] = oscillateurAmorti(f,duration,f0,alpha,fs)
    k = 1/fs;
    Nt = floor(duration*fs);
    if length(f) > Nt
        error("Force d'excitation plus longue que la durée")
    end
    t = (0:Nt-1)*k;
    f = [f zeros(1,Nt - length(f))]; % zero-padding
    u = zeros(size(t));
    omega0 = 2*pi*f0;
    
    % paramètres dérives
    % (Calcul préalable à la boucle pour optimiser)
    c1 = 2-k^2*omega0^2;
    c2 = (alpha*k-1);
    c3 = k^2;
    c4 = 1 + alpha*k;
    
    % boucle temporelle
    for n = 2 : Nt-1
        u(n+1) = (u(n)*c1 + u(n-1)*c2 + f(n)*c3)/c4;
    end
end


