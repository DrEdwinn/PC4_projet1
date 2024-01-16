function [l_harmo] = readHarmonique(partition)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
l_mesure = partition.part.measure;
l_harmo = [];
for m = 1:length(l_mesure)
    l_note = l_mesure(m).note;
    for n = 1:length(l_note)
        if ~ismissing(l_note(n).notations)
            if isfield(l_note(n).notations.technical,"harmonic")
                l_harmo = [l_harmo "harmonic"];
            else
                l_harmo = [l_harmo "rien"];
            end
        else
            l_harmo = [l_harmo "rien"];
        end
    end
end
end