function [l_dur] = readDurations(partition)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
l_mesure = partition.part.measure;
l_dur = [];
for m = 1:length(l_mesure)
    l_note = l_mesure(m).note;
    for n = 1:length(l_note)
        l_dur = [l_dur l_note(n).duration];
    end
end
end