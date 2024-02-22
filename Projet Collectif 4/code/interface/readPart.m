function [excitation, string] = readPart(time, part)
% READPART - Lit les paramètres d'excitation et de corde depuis une portée
% 
% Cette fonction renvoie deux structures de données correspondant aux
% paramètres physiques de l'excitation et de la corde lue depuis une
% matrice contenant les informations écrites dans la portée d'une
% partition.
% 
% Syntaxe
%   [excitation, string] = READPART(time, part, preset)
% 
%   excitation : structure de données contenant les valeurs des paramètres
%       physiques de l'excitation au cours du temps
% 
%   string : structure de données contenant les valeurs des paramètres
%       physiques de la corde au cours du temps
% 
%   time : structure de données contenant les paramètres généraux liés au
%       temps
% 
%   part : matrice contenant les informations d'une portée de la partition
% 
%   preset : structure de données contenant les paramètres physiques de
%       départ

% paramètres de la corde
presetName = part{1, index.PRESET_NAME};
if presetName == ""
    presetName = 'default';
end
preset = loadPreset(presetName);
string.L = preset.L;
string.sigma0 = preset.sigma0;
string.sigma1 = preset.sigma1;
string.xrecord = preset.xrecord;
string.sigmar = preset.sigmar;
string.xratt = preset.xratt;
string.omegaratt = preset.omegaratt;
string.mratt = preset.mratt;

% paramètres d'excitation
excitation.xe = preset.xe;
excitation.fPluck = zeros(size(time.t));
excitation.fBow = zeros(size(time.t));
excitation.amp = preset.amp;
excitation.type = preset.exType;
excitation.pluckDuration = preset.pluckDuration;
excitation.bowForce = preset.bowForce;
excitation.bowVelocity = preset.bowVelocity * ones(size(time.t));
excitation.bowAttack = preset.bowAttack;
excitation.bowModFreq = preset.bowModFreq;
excitation.bowModAmount = preset.bowModAmount;

% calcul de la fonction d'excitation à partir des notes
[excitation.fPluck, excitation.fBow] = getExFunction(time,part,excitation); % fonction à renommer peut être

% calcul de la position du ressort à partir des notes
f0 = preset.gamma/string.L/2; 
freqRatio = f0 * getNotesFrequency(time, part).^-1;
if isfield(preset,'xp')
    string.xp = preset.xp * ones(size(time.t)); 
else
    string.xp = min(1,freqRatio) * string.L;
end

% calcul des fonctions des paramètres physique au cours du temps
if excitation.type=='mouette'
    string.gamma = gammaMouette(time, part, preset.gamma);
else
    detune = getParameter(time, part, index.DETUNE_VALUE, 0);
    string.gamma = preset.gamma * 2.^(detune/12);
end
string.kappa = getParameter(time, part, index.KAPPA_VALUE, preset.kappa);
string.xr = getParameter(time, part, index.XR_VALUE, preset.xr);
string.omega0 = getParameter(time, part, index.OMR0_VALUE, preset.omega0);
string.omega1 = getParameter(time, part, index.OMR1_VALUE, preset.omega1);
end

