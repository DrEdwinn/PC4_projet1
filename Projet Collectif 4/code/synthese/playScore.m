function [out] = playScore(score, fs)
%PLAYSCORE - Renvoie l’audio synthétisé par les modèles de cordes à partir d’une partition.
%
% out : une matrice de vecteurs audio où chaque ligne correspond à la synthèse d’une portée de la partition.
% 
% score : structure de donnée contenant les information de la partition 
% 
% fs : fréquence d’échantillonage


% paramètres temporels
time.fs = fs;
time.t = (0:score.duration*time.fs-1)/time.fs;

% pré-allocation
out = zeros(score.partCount,length(time.t));


% parcours des portées de la partition
for iPart = 1:score.partCount
    part = getPart(score.notes,iPart);              % extraction de la i-ème portée
    [excitation, string] = readPart(time, part);    % conversion des informations de la portée
                                                    % en paramètres physiques
    string.xp = zeros(size(time.t));
    out(iPart,:) = (stiffString(time, excitation, string)); % synthèse
end
end

