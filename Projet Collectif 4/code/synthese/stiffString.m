function [out] = stiffString(T, ex, s)
%STIFFSTRING - Renvoie le son synthétisé du modèle de la corde préparée à partir des paramètres physiques
% 
%  out : vecteur audio contenant le son synthétisé
% 
%  T : structure de données contenant les paramètres généraux liés au
%      temps
% 
%  s : structure de données contenant les valeurs des paramètres
%      physiques de la corde au cours du temps
% 
%  ex : structure de données contenant les valeurs des paramètres
%      physiques de l'excitation

    flagPlot = 0; % plot la corde si = 1
    k = 1/T.fs;
    Nt = length(T.t);

    
    
    gamma = s.gamma;
    kappa = s.kappa;
    L = s.L;
    sigma0 = s.sigma0;
    sigma1 = s.sigma1;

    % condition de stabilité
    h = sqrt((max(gamma).^2 * k.^2 + sqrt(max(gamma).^4 *k.^4 + 16 * max(kappa).^2 *k.^2))/2)*1.1;
    Nx = floor(L/h);
    h = L/Nx;
    x = (0:Nx)*h;
    
    
    u = zeros(3,length(x));     % déplacement de la corde, les lignes 1 à 3 
                                % correspondant aux temps n-2 à n

   
    u_ratt = zeros(size(T.t));  %déplacement du rattling element

    % paramètre ressorts
    omegap_0 = 25000;
    omegap_1 = 0;
    sigma_p = 0;
    omegar_0 = s.omega0;
    omegar_1 = s.omega1;
    sigma_r = s.sigmar;
    
    % excitation
    x_e = ex.xe;
    fpluck = ex.fPluck;
    fbow = zeros(size(T.t));

    % paramètres archet
    alpha_bow = 100;
    vb = ex.bowVelocity;
    
    % paramètres rattling element
    x_ratt = s.xratt;            % position
    omega_ratt = s.omegaratt;       % raideur (non-linéaire)
    alphaR = 3;             % exposant de la raideur
    epsilon_ratt = 5e-8;  % longueur
    M = s.mratt;                  % rapport masse de la corde 
                            % sur masse du rattling element

    

    % paramètres dérivés
    lambda = k*gamma/h;
    eps = 10^(-6);
    cb = sqrt(2*alpha_bow) * exp(1/2);
    C1 = 2 - 2*lambda.^2 - 6 * k^2*kappa.^2/h^4 - 4 * sigma1 *k/h^2;
    C2 = -1 + sigma0 *k + 4 * sigma1 * k /h^2;
    C3 = lambda.^2 + 4 * k^2 * kappa.^2/ h^4 + 2 * sigma1 *k/h^2;
    C4 = kappa.^2*k^2 /h^4;
    C5 = 2 * sigma1 * k /h^2;

    out = zeros(size(T.t));
    
    
    %boucle temporelle
    for n = 2:Nt-1
        %%%%%%%%%%%%%%% operateurs d'interpolation %%%%%%%%%%%%%%%%
        %%% position ressort doigts 
        x_p = s.xp(n);
        l_p = floor(x_p/h)+1;           %partie entière
        alpha_p = x_p/h - (l_p-1);      %partie fractionnaire
        l_p = min(l_p,Nx-2);
        l_p = max(l_p,3);
        J_l1_p = zeros(size(x));
        J_l1_p(l_p) = (1-alpha_p)/h;
        J_l1_p(l_p+1) = alpha_p/h;

        %%% position ressort perturbateur
        x_r = s.xr(n);
        l_r = floor(x_r/h)+1;           %partie entière
        alpha_r = x_r/h - (l_r-1);      %partie fractionnaire
        l_r = min(l_r,Nx-2);
        l_r = max(l_r,3);
        J_l1_r = zeros(size(x));
        J_l1_r(l_r) = (1-alpha_r)/h;
        J_l1_r(l_r+1) = alpha_r/h;        
        
        %%% position d'excitation
        l_e = floor(x_e/h)+1;           %partie entière
        alpha_e = x_e/h - (l_e-1);      %partie fractionnaire
        l_e = min(l_e,Nx-2);
        l_e = max(l_e,3);
        J_l1 = zeros(size(x));
        J_l1(l_e) = (1-alpha_e)/h;
        J_l1(l_e+1) = alpha_e/h;

        %%% position rattling element
        l_ratt = floor(x_ratt/h)+1;         %partie entière
        alpha_ratt = x_ratt/h - (l_ratt-1); %partie fractionnaire
        J_l1_ratt = zeros(size(x));
        J_l1_ratt(l_ratt) = (1-alpha_ratt)/h;
        J_l1_ratt(l_ratt+1) = alpha_ratt/h;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Calcul force ressort doigts

        eta_pas = (1- alpha_p) * u(1,l_p) + alpha_p * u(1,l_p+1);
        eta_pre = (1- alpha_p) * u(2,l_p) + alpha_p * u(2,l_p+1);
        eta_pas_lp1 = (1- alpha_p) * u(1,l_p+1) + alpha_p * u(1,l_p+2);
        eta_pre_lp1 = (1- alpha_p) * u(2,l_p+1) + alpha_p * u(2,l_p+2);
        eta_pas_lm1 = (1- alpha_p) * u(1,l_p-1) + alpha_p * u(1,l_p);
        eta_pre_lm1 = (1- alpha_p) * u(2,l_p-1) + alpha_p * u(2,l_p);
        eta_pre_lp2 = (1- alpha_p) * u(2,l_p+2) + alpha_p * u(2,l_p+3);
        eta_pre_lm2 = (1- alpha_p) * u(2,l_p-2) + alpha_p * u(2,l_p-1);

        %paramètres dérivés
        a = 1 + k*sigma0 + k^2/h*(omegap_0^2/2 + omegap_1^4/2 + sigma_p/k);
        b = -1 + k*sigma0 -2*sigma1*k*(eta_pas_lp1 -2*eta_pas + eta_pas_lm1)...
            -k^2/h*(omegap_0^2/2 + omegap_1^4/2 - sigma_p/k);
        c = eta_pre*(2 + 2*sigma1*k*(eta_pre_lp1 -2*eta_pre + eta_pre_lm1)...
            + k^2*gamma(n)^2*(eta_pre_lp1 -2*eta_pre + eta_pre_lm1)...
            - k^2*kappa(n)^2*(eta_pre_lp2 -4*eta_pre_lp1 + 6*eta_pre -4*eta_pre_lm1 + eta_pre_lm2));

        eta_fut = c/a + b/a*eta_pas;
        Fp = -omegap_0^2/2*(eta_fut+eta_pas) -omegap_1^4/2*eta_pre^2*(eta_fut+eta_pas) -sigma_p/k*(eta_fut - eta_pas);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Calcul force ressort perturbateur

        eta_pas = (1- alpha_r) * u(1,l_r) + alpha_r * u(1,l_r+1);
        eta_pre = (1- alpha_r) * u(2,l_r) + alpha_r * u(2,l_r+1);
        eta_pas_lp1 = (1- alpha_r) * u(1,l_r+1) + alpha_r * u(1,l_r+2);
        eta_pre_lp1 = (1- alpha_r) * u(2,l_r+1) + alpha_r * u(2,l_r+2);
        eta_pas_lm1 = (1- alpha_r) * u(1,l_r-1) + alpha_r * u(1,l_r);
        eta_pre_lm1 = (1- alpha_r) * u(2,l_r-1) + alpha_r * u(2,l_r);
        eta_pre_lp2 = (1- alpha_r) * u(2,l_r+2) + alpha_r * u(2,l_r+3);
        eta_pre_lm2 = (1- alpha_r) * u(2,l_r-2) + alpha_r * u(2,l_r-1);

        %paramètres dérivés
        a = 1 + k*sigma0 + k^2/h*(omegar_0(n)^2/2 + omegar_1(n)^4/2 + sigma_r/k);
        b = -1 + k*sigma0 -2*sigma1*k*(eta_pas_lp1 -2*eta_pas + eta_pas_lm1)...
            -k^2/h*(omegar_0(n)^2/2 + omegar_1(n)^4/2 - sigma_r/k);
        c = eta_pre*(2 + 2*sigma1*k*(eta_pre_lp1 -2*eta_pre + eta_pre_lm1)...
            + k^2*gamma(n)^2*(eta_pre_lp1 -2*eta_pre + eta_pre_lm1)...
            - k^2*kappa(n)^2*(eta_pre_lp2 -4*eta_pre_lp1 + 6*eta_pre -4*eta_pre_lm1 + eta_pre_lm2));

        eta_fut = c/a + b/a*eta_pas;
        Fr = -omegar_0(n)^2/2*(eta_fut+eta_pas) -omegar_1(n)^4/2*eta_pre^2*(eta_fut+eta_pas) -sigma_r/k*(eta_fut - eta_pas);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%calcul de la force de l'archet

        if strcmp(ex.type,'bow') || strcmp(ex.type,'mouette') || strcmp(ex.type,'bowAndPluck')
            eta_pas = (1- alpha_e) * u(1,l_e) + alpha_e * u(1,l_e+1);
            eta_pre = (1- alpha_e) * u(2,l_e) + alpha_e * u(2,l_e+1);
            eta_pas_lp1 = (1- alpha_e) * u(1,l_e+1) + alpha_e * u(1,l_e+2);
            eta_pre_lp1 = (1- alpha_e) * u(2,l_e+1) + alpha_e * u(2,l_e+2);
            eta_pas_lm1 = (1- alpha_e) * u(1,max(l_e-1,1)) + alpha_e * u(1,l_e);
            eta_pre_lm1 = (1- alpha_e) * u(2,max(l_e-1,1)) + alpha_e * u(2,l_e);
            eta_pre_lp2 = (1- alpha_e) * u(2,l_e+2) + alpha_e * u(2,l_e+3);
            eta_pre_lm2 = (1- alpha_e) * u(2,max(l_e-2,1)) + alpha_e * u(2,max(l_e-1,1));
            
            % constante contenant les éléments connus pour l'algo de N-R
            b = (2*k * vb(n) + eta_pas) * (1 + sigma0 * k ) - eta_pre * C1(n) ...
                - eta_pas * C2 - (eta_pre_lp1 + eta_pre_lm1) * C3(n) ...
                + C4(n) * (eta_pre_lp2+eta_pre_lm2) + C5 * (eta_pas_lp1 + eta_pas_lm1);
            n_iter = 1;         %compteur pour l'algo de N-R
            dif = 1;            %initialisation de la différence 

            %algorithme de Newton-Raphson pour le calcul de la vitesse relative
            %entre l'archet et la corde
            v_rel = zeros(1,101);
            F = zeros(1,100);
            while dif>eps
                F(n_iter) = k^2/h * ex.fBow(n) * cb * v_rel(n_iter) * exp(-alpha_bow*v_rel(n_iter)^2) ...
                    + 2*k*v_rel(n_iter) * (1 + sigma0 * k) + b;
                dF = k^2/h * ex.fBow(n) *cb * exp(-alpha_bow*v_rel(n_iter)^2) * (1 - 2*alpha_bow*v_rel(n_iter)^2) ...
                    + 2*k * (1 + sigma0 * k);
                v_rel_new = v_rel(n_iter) - F(n_iter)/dF;
                dif = abs(v_rel_new-v_rel);
                n_iter = n_iter+1;
                v_rel(n_iter) = v_rel_new;
                if n_iter>100
                    figure(1)
                    plot(v_rel(1:100),F,'.');
                    error('nombre d iteration trop important');
                end
            end

            %calcul de la force de l'archet sur la corde
            fbow(n) = - ex.fBow(n) * ex.bowForce * cb * v_rel_new * exp(-alpha_bow*v_rel_new^2);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% calcul force rattling element

        eta_pre = (1- alpha_ratt) * u(2,l_ratt) + alpha_ratt * u(2,l_ratt+1);

        if eta_pre - u_ratt(n) >= epsilon_ratt/2
            F_ratt = -omega_ratt^(alphaR+1)*(eta_pre - u_ratt(n) + epsilon_ratt/2)^(alphaR+1);
        elseif abs(eta_pre - u_ratt(n)) < epsilon_ratt/2
            F_ratt = 0;
        else
            F_ratt = omega_ratt^(alphaR+1)*(eta_pre - u_ratt(n) + epsilon_ratt/2)^(alphaR+1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Calcul déplacement de la corde

        % Conditions limite supportées
        l=2;
        u(3,l) = (u(2,l) * C1(n) ...
            + u(1,l)*C2 ...
            + u(2,l+1)*C3(n) ...
            - (u(2,l+2)-u(2,l)) * C4(n) ...
            - (u(1,l+1)) * C5 ...
            + J_l1(l)*(fpluck(n) + fbow(n))*ex.amp*k^2 + J_l1_p(l)*Fp*k^2 ...
            + J_l1_r(l)*Fr*k^2 + J_l1_ratt(l)*F_ratt*k^2) ...
            / (1 + sigma0 * k);
        l=Nx;
        u(3,l) = (u(2,l) * C1(n) ...
            + u(1,l)*C2 ...
            + u(2,l-1)*C3(n) ...
            - (u(2,l-2)-u(2,l)) * C4(n) ...
            - (u(1,l-1)) * C5 ...
            + J_l1(l)*(fpluck(n) + fbow(n))*ex.amp*k^2 + J_l1_p(l)*Fp*k^2 ...
            + J_l1_r(l)*Fr*k^2 + J_l1_ratt(l)*F_ratt*k^2) ...
            / (1 + sigma0 * k);

        % Parcours des échantillons de la cordes        
        u(3,3:Nx-1) = (u(2,3:Nx-1) * C1(n) ...
            + u(1,3:Nx-1) * C2 ...
            + (u(2,4:Nx) + u(2,2:Nx-2)) * C3(n) ...
            - (u(2,5:Nx+1)+u(2,1:Nx-3)) * C4(n) ...
            - (u(1,4:Nx)+u(1,2:Nx-2)) * C5 ...
            + J_l1(3:Nx-1)*(fpluck(n) + fbow(n))*ex.amp*k^2 + J_l1_p(3:Nx-1)*Fp*k^2 ...
            + J_l1_r(3:Nx-1)*Fr*k^2 + J_l1_ratt(3:Nx-1)*F_ratt*k^2) ...
            / (1 + sigma0 * k);        
        
        %sortie sur un point d'écoute
        out(n) = u(2,floor(s.xrecord/h)+2);
        
        if flagPlot ==1
            figure(1)
            %pause(.05)
            plot(x,u(3,:),'b')
            hold on
            plot(x_ratt,u_ratt(n),'or')
            hold off
            ylim([-3*10^(-10) 3*10^(-10)])
            drawnow
        end
        %mise à jour des déplacements itération suivante
        u_ratt(n+1) = 2*u_ratt(n) -u_ratt(n-1) -k^2*M*F_ratt;
        u(1,:) = u(2,:);
        u(2,:) = u(3,:);
    end
end



