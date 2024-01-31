function [parameters] = randParameters()
%RANDPARAMETERS Summary of this function goes here
%   Detailed explanation goes here

% paramètres de corde par défault
parameters.L = .005 * 10000^rand() + 0.995;
parameters.gamma = rand() * (5000 - 10) + 10;
parameters.kappa =  7*100^rand() - 1;
parameters.sigma0 = .5 * 100^rand()-1;
parameters.sigma1 = rand() * (0.1 - 0.00001) + 0.00001;
parameters.xrecord = 0.1;
parameters.xp = .9 * rand() + 0.1 ;
parameters.omega0 = 1000000*rand()-1;
parameters.omega1 = 1000000*rand()-1;

% paramètres d'excitation par défault
parameters.xe = rand();
parameters.amp = rand() * (100000 - 0.1) + 0.1;
parameters.pluckDuration = rand() * (0.01 - 0.0001) + 0.0001;
parameters.bowVelocity = rand() * (0.3 - 0.03) + 0.03;
parameters.bowForce = rand() * 149 + 1;

exTypeList = ["pluck" "bow"];
parameters.exType = exTypeList(floor(rand()*length(exTypeList))+1);

end

