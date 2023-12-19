% matlab script idealbarfd.m
% finite difference scheme for the ideal bar equation 
% clamped/pivoting boundary conditions
% raised cosine initial conditions

%%%%%% begin global parameters

SR = 44100;                 % sample rate (Hz)
K = 10;                     % stiffness parameter
TF = 1;                     % duration of simulation (s)
ctr = 0.7; wid = 0.1;       % center location/width of excitation
u0 = 1; v0 = 0;             % maximum initial displacement/velocity
mu = 0.5;                   % scheme free parameter
rp = 0.85;                  % position of readout (0-1)
bc = [2 2];                 % boundary condition type, [left right] with
                            % 1: clamped, 2: pivoting

%%%%%% end global parameters

% begin derived parameters

k = 1/SR;                   % time step
NF = floor(SR*TF);          % duration of simulation (samples)

% stability condition/scheme parameters

h = sqrt(K*k/mu); N = floor(1/h); h = 1/N; mu = K*k/h^2; 
s0 = 2*(1-3*mu^2); s1 = 4*mu^2; s2 = -mu^2;

% readout interpolation parameters

rp_int = 1+floor(N*rp);     % rounded grid index for readout 
rp_frac = 1+rp/h-rp_int;    % fractional part of readout location

% create raised cosine

xax = [0:N]'*h;
ind = sign(max(-(xax-ctr-wid/2).*(xax-ctr+wid/2),0)); 
rc = 0.5*ind.*(1+cos(2*pi*(xax-ctr)/wid));

% initialize grid functions and output

u2 = u0*rc; u1 = (u0+k*v0)*rc; u = zeros(N+1,1); out = zeros(NF,1); 

%%%%%% start main loop
for n=3:NF
    % scheme calculation (interior)
    u(3:N-1) = -u2(3:N-1)+s0*u1(3:N-1)+s1*(u1(2:N-2)+u1(4:N))...
        +s2*(u1(1:N-3)+u1(5:N+1));
    % calculations at boundary points 
    if(bc(1)==2)
        u(2) = -u2(2)+(s0-s2)*u1(2)+s1*u1(3)+s2*u1(4);
    end
    if(bc(2)==2)
       u(N) = -u2(N)+(s0-s2)*u1(N)+s1*u1(N-1)+s2*u1(N-2);
    end
    out(n) = (1-rp_frac)*u(rp_int)+rp_frac*u(rp_int+1); % readout
    u2 = u1; u1 = u; % update

end

%%%%% end main loop

% plot output waveform

plot([0:NF-1]*k, out, 'k');
xlabel('t'); ylabel('u'); title('Ideal Bar Equation: FD Output'); 
axis tight

% play sound

soundsc(out,SR);

