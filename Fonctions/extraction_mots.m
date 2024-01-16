function [l_word] = extraction_mots(partition)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
l_mesure = partition.part.measure;
l_word = [];
l_dynamic = [];
for m = 1:length(l_mesure)
    l_direction = l_mesure(m).direction;
    for d = 1:length(l_direction)
        if isfield(l_direction(d).direction_type,"words")
            l_word = [l_word l_direction(d).direction_type.words];
        end
    end
end
end