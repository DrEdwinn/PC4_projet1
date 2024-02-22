function [] = savePreset(parameters, name)
%SAVEPRESET - Écrit une structure de paramètres dans un fichier
%
% Cette fonction crée fichier à partir d'une structure de paramètres
%
% Syntaxe
%   SAVEPRESET(parameters, name)
% 
%   parameters : structure de données correspondant aux valeurs de
%       paramètres par défault de la corde et de l'excitation.
%
%   name : chaine de caractère correspondant à un nom de fichier
 

name = erase(name, ' '); % suppression des espaces dans le nom du fichier
filename = "../../presets/" + name + ".xml"; % ajout de l'extension .xml
writestruct(parameters,filename);
end

