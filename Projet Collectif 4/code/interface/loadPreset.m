function [parameters] = loadPreset(filename)
% LOADPRESET - Lit un jeu de paramètres depuis un fichier
% 
% Cette fonction crée une structure de paramètres lue depuis le fichier
% spécifié.
% 
% Syntaxe
%   parameters = LOADPRESET(filename)
% 
%   filename : chaine de caractère correspondant à un nom de fichier
% 
%   parameters : structure de données correspondant au jeu de
%       paramètres par défault de la corde et de l'excitation.

if ~contains(filename, ".xml")
    filename = filename + ".xml"; % ajout de l'extension .xml si absente
end
filepath = char(which(filename));
parameters = readstruct(filepath);
end

