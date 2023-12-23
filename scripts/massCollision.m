% Modèle : eq (11.a) https://arxiv.org/pdf/1405.2589.pdf

clc
fs = 48000;
k = 1/fs;
duration = 3;
Nt = floor(fs*duration);
t = (0:Nt-1)*k;

M = 1; % massez moi les épaules
K = 2000;
alpha = 1; % ordre de la loi de PUISSANCE
yc = 200; % position de la collision

u = zeros(size(t));
f = zeros(size(t));

u(1) = 0;
u(2) = .01;
udelta = u - yc;

eps = 10^(-6);

for n=3:Nt-1
    a = u(n-2);
    m = k^2/M;
    r = yc; % estimation de q0 pour newton-raphson
    b = -2*u(n-1) + 2*u(n-2);
    n_iter = 0;
    delta = 1;
        
    while delta>eps % Newton-fucking-Raphson
        % if r<10^(-12)
        %     G = r + m*K*(a - yc)^alpha + b;
        %     dG = 1 - m*alpha*K*(a -yc)^(alpha-1);
        % else
        G = r + m/r*(phi(r + a,yc,K,alpha) - phi(a,yc,K,alpha)) + b;
        dG = 1 - m/r^2*(phi(r+a,yc,K,alpha) - phi(a,yc,K,alpha)) + m/r*K*(r+a-yc)^alpha;
        %end
        n_iter = n_iter +1;
        r_new = r - G/dG;
        delta = abs(r_new - r);
        r = r_new;
        if n_iter > 100
            error("trop d'iterations");
        end
    end
    u(n) = r + a;
    udelta(n) = u(n) - yc;
    f(n-1) = 1/2/k*(phi(u(n),yc,K,alpha) + phi(u(n-2),yc,K,alpha));
end

figure(1);
yyaxis left;
plot(t,u);
ylabel("déplacement")
yyaxis right
plot(t,f);
ylabel("force")
xlabel("temps (s)")

function [w] = phi(u, yc, K, alpha)
    v = u - yc;
    if v <= 0 % partie positive
        v = 0;
    end
    w = K*v^(alpha+1)/(alpha+1);
end