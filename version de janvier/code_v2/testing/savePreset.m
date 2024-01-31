function [] = savePreset(parameters, name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
name = erase(name, ' ');
filename = "../../presets/" + name + ".xml";
writestruct(parameters,filename);
end

