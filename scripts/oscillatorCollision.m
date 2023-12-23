% Modèle : eq (11.a) https://arxiv.org/pdf/1405.2589.pdf

clc
fs = 48000;
k = 1/fs;
duration = 3;
Nt = floor(fs*duration);
t = (0:Nt-1)*k;

M = 1; % massez moi les épaules
K = 100000000;
alpha = 1; % ordre de la loi de PUISSANCE
yc = .1; % position de la collision
sigma0 = 20; % coefficient d'amortissement
f0 = 60;
omega0 = 2*pi*f0;

u = zeros(size(t));
u(1) = -10;
u(2) = -9;
fc = zeros(size(t)); % force d'interaction
fex = zeros(size(t)); % force d'excitation
%fex(1:100) = hanning(100);

eps = 10^(-6);

for n=3:Nt-1
    a = u(n-2);
    m = k^2/M/(1+sigma0*k);
    r = 0; % estimation de q0 pour newton-raphson
    b = (u(n-1)*(k^2*omega0^2 - 2) + 2*u(n-2))/(1 + sigma0*k) - m*fex(n-1);
    n_iter = 0;
    delta = 1;
        
    while delta>eps % Newton-fucking-Raphson
        if abs(r)<10^(-12)
            G = r + m*K*positivePart(a - yc)^alpha + b;
            dG = 1 - m*alpha*K*positivePart(a -yc)^(alpha-1);
        else
            G = r + m/r*(phi(r + a,yc,K,alpha) - phi(a,yc,K,alpha)) + b;
            dG = 1 - m/r^2*(phi(r+a,yc,K,alpha) - phi(a,yc,K,alpha)) + m/r*K*positivePart(r+a-yc)^alpha;
        end
        n_iter = n_iter +1;
        r_new = r - G/dG;
        delta = abs(r_new - r);
        r = r_new;
        if n_iter > 100
            error("trop d'iterations");
        end
    end
    u(n) = r + a;
    fc(n-1) = 1/2/k*(phi(u(n),yc,K,alpha) + phi(u(n-2),yc,K,alpha));
end

figure(1);
yyaxis left;
plot(t,u);
ylabel("déplacement")
%yyaxis right
%plot(t,fc);
%ylabel("force")
xlabel("temps (s)")
xlim([0 0.1])
sgtitle("Collision d'un oscillateur à une paroi rigide")
%soundsc(u,fs);
function [w] = phi(u, yc, K, alpha)
    w = K*positivePart(u - yc)^(alpha+1)/(alpha+1);
end

