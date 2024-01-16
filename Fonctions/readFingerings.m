function [l_pos] = readFingerings(partition)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
l_mesure = partition.part.measure;
l_pos = [];
for m = 1:length(l_mesure)
    l_note = l_mesure(m).note;
    for n = 1:length(l_note)
        if ~ismissing(l_note(n).notations)
            if isfield(l_note(n).notations.technical, "fingering")
                l_pos = [l_pos l_note(n).notations.technical.fingering];
            else
                l_pos = [l_pos 0];
            end
        else
            l_pos = [l_pos 0];
        end
    end
end
end