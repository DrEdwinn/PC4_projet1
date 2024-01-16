function [l_notehead] = readNoteheads(partition)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
l_mesure = partition.part.measure;
l_notehead = [];
for m = 1:length(l_mesure)
    l_note = l_mesure(m).note;
    for n = 1:length(l_note)
        if l_note(n).notehead == "x"
            l_notehead = [l_notehead "frapper"];
        elseif l_note(n).notehead == "triangle"
                l_notehead = [l_notehead "pincer"];
        elseif l_note(n).notehead == "slashed"
                l_notehead = [l_notehead "frotter"];
        else
            l_notehead = [l_notehead "silence"];
        end
    end
end
end