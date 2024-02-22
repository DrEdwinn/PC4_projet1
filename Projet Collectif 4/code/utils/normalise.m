function [y] = normalise(x, part, time)
%NORMALISE - Renvoie un vecteur audio où le volume de chaque note est normalisée.
% 
% Cette fonction fait face au problème où de grandes différences d’amplitude était observée selon la hauteur de la note, ce qui en rendait certaines inaudibles.
% 
%  y : vecteur audio de sortie normalisée
% 
%  x : vecteur audio d’entrée
% 
%  time : structure de données contenant les paramètres généraux liés au
%      temps
% 
%  part : matrice contenant les informations d'une portée de la partition

y = zeros(size(x));
noteOnset = [part{:,3}];
noteDuration = [part{:,4}];
isRest = [part{:,5}]==0;
L = length(noteOnset);
for i = 1:L
    if ~isRest(i) % si la note n'est pas un silence
        onset = floor(noteOnset(i)*time.fs)+1;
        offset = floor(noteDuration(i)*time.fs);
        x_i = x(onset:onset+offset-1);
        if any(x_i) % 
            y(onset:onset+offset-1) = x_i/max(abs(x_i));
        end
    end
end