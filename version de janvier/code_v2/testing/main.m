clc
% import de la partition
filename = input("Entrer le nom du fichier xml : ",'s');
[Score, duration, partCount] = parseMusicXml(filename);

% paramètres généraux
time.fs = 48000; 
time.t = (0:duration*time.fs-1)/time.fs;
out = zeros(partCount,length(time.t));

% importation du preset
parameters = loadPreset("basic_string");
parameters.gamma = 330;

% synthèse;
for iPart = 1:partCount
    partInd = partIndex(Score,iPart);
    part = Score(partInd);
    part = reshape(part,length(part)/width(Score),width(Score));
    [excitation, string] = readPart(time, part, parameters);
    out(iPart,:) = stiffString(time, excitation, string);
end

voiceVol = ones(1,partCount);
if ~any(isnan(out))
    soundsc(sum(voiceVol*out,1) ,time.fs)
    if any(out(1,:))
        out(1,:) = out(1,:)/max(abs(out(1,:)));
    else
        warning("PAS DE SON");
    end
else
    warning('INPUT CONTAINS NAN VALUE WHICH CAN RESULT IN UNEXPECTED BEHAVIOR')
end

clear("parameters.xp");














