function [parameters] = loadPreset(filename)
%LOADPRESET Summary of this function goes here
%   Detailed explanation goes here

if ~contains(filename, ".xml")
    filename = filename + ".xml";
end
filepath = "../../presets/";
parameters = readstruct(filepath + filename);
end

