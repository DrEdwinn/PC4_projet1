function [hertz] = readFrequencies(partition)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
l_mesure = partition.part.measure;
l_step = [];
l_oct = [];
for m = 1:length(l_mesure)
    l_note = l_mesure(m).note;
    for n = 1:length(l_note)
        if ~ismissing(l_note(n).pitch)
            l_step = [l_step  l_note(n).pitch.step];
            l_oct = [l_oct l_note(n).pitch.octave];
        else
            l_step = [l_step "0"];
            l_oct = [l_oct 0];
        end
    end
end

hertz = [];
for i = 1:length(l_oct)
    if l_step(i) == "C"
        hertz = [hertz 32.70*2^(l_oct(i))];
    end
    if l_step(i) == "D"
        hertz = [hertz 36.71*2^(l_oct(i))];
    end
    if l_step(i) == "E"
        hertz = [hertz 41.2*2^(l_oct(i))];
    end
    if l_step(i) == "F"
        hertz = [hertz 43.65*2^(l_oct(i))];
    end
    if l_step(i) == "G"
        hertz = [hertz 49*2^(l_oct(i))];
    end
    if l_step(i) == "A"
        hertz = [hertz 55*2^(l_oct(i))];
    end
    if l_step(i) == "B"
        hertz = [hertz 61.74*2^(l_oct(i))];
    end
    if l_step(i) == "0"
        hertz = [hertz 0];
end
end