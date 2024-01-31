function [out] = singleNote(time, parameters)
%SINGLENOTE Summary of this function goes here
%   Detailed explanation goes here

% import de la partition
[Score, duration] = parseMusicXml('singlenote.musicxml');
partInd = partIndex(Score,1);
part = Score(partInd);
part = reshape(part,length(part)/width(Score),width(Score));

% paramètres généraux
time.t = (0:duration*time.fs-1)/time.fs;

%synthèse
[excitation, string] = readPart(time, part, parameters);
out = stiffString(time, excitation, string);
end

