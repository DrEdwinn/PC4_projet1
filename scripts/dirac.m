function [f] = dirac(t0,fs)
    Nt = floor(t0*fs);
    f = zeros(1,Nt);
    f(Nt) = 1;
end