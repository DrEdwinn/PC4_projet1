% matlab script sho.m
% finite difference scheme for simple harmonic oscillator
%%%%%% begin global parameters
SR = 44100;
f0 = 1000;
TF = 1.0;
u0 = 0.3;
v0 = 0.0;
%%%%%% end global parameters
% check that stability condition is
% sample rate (Hz)
% fundamental frequency (Hz) % duration of simulation (s) % initial displacement
% initial velocity
if(SR<=pi*f0)
    error('Stability condition violated');
end
% derived parameters
k = 1/SR;
coef = 2-k^2*(2*pi*f0)^2;
NF = floor(TF*SR);
% initialize state of scheme
u1 = u0+k*v0;
u2 = u0;
% time step
% scheme update coefficient
% duration of simulation (samples)
% last value of time series
% one before last value of time series
% initialize readout
out = zeros(NF,1); out(1) = u2; out(2) = u1; %%%%%% start main loop
for n=3:NF
    u=coef*u1-u2;
out(n) = u; u2=u1;u1=u;
end
%%%%%% end main loop
% play sound
soundsc(out,SR);
% plot output waveform
% difference scheme calculation
% read value to output vector
% update state of difference scheme
plot([0:NF-1]*k, out, 'k');
xlabel('t'); ylabel('u'); title('SHO: Scheme Output'); axis tight