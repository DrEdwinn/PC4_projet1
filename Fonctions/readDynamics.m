function [l_dyn] = readDynamics(partition)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
l_mesure = partition.part.measure;
l_dyn = [];
for m = 1:length(l_mesure)
    l_direction = l_mesure(m).direction;
    for d = 1:length(l_direction)
        if isfield(l_direction(d).direction_type,"dynamics")
            if isfield(l_direction(d).direction_type.dynamics,"f")
                l_dyn = [l_dyn "f"];
            elseif isfield(l_direction(d).direction_type.dynamics,"mp")
                l_dyn = [l_dyn "mp"];
            end
        end
    end
end
end