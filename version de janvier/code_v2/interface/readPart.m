function [excitation, string] = readPart(time, part, parameters)
%READPART Summary of this function goes here
%   Detailed explanation goes here

%string parameters
string.L = parameters.L;
string.sigma0 = parameters.sigma0;
string.sigma1 = parameters.sigma1;
string.xrecord = parameters.xrecord;
string.sigmar = 0;

%excitation parameters
excitation.amp = parameters.amp;
excitation.type = parameters.exType;
excitation.pluckDuration = parameters.pluckDuration;
excitation.bowForce = parameters.bowForce;
excitation.bowVelocity = parameters.bowVelocity * ones(size(time.t));

%compute excitation function based on notes
excitation.f = getExFunction(time, part, excitation);

%compute spring position based on notes
f0 = parameters.gamma/string.L/2; 
freqRatio = f0 * getNotesFrequency(time, part).^-1;
if isfield(parameters,'xp')
    string.xp = parameters.xp * ones(size(time.t)); 
else
    string.xp = min(1,freqRatio) * string.L;
end

%compute physical parameters time-vectors
string.gamma = getParameter(time, part, 6, parameters.gamma);
string.kappa = getParameter(time, part, 8, parameters.kappa);
excitation.xe = getParameter(time, part, 10, parameters.xe);
string.xr = getParameter(time, part, 12, parameters.xr);
string.omega0 = getParameter(time, part, 14, parameters.omega0);
string.omega1 = getParameter(time, part, 16, parameters.omega1);
end

