function [freq] = getNotesFrequency(time, Part)
%GETNOTESFREQUENCY returns a vector of same size as time vector t which 
%contains the frequencies written in thes part over time.
%   à compléter...
currFreq = 0;
freq = zeros(size(time.t));
noteOnset = Part(:,3);
noteDuration = Part(:,4);
notePitch = Part(:,5);
L = length(noteOnset);
for i=1:L
    onset = floor(noteOnset(i)*time.fs)+1;
    offset = floor(noteDuration(i)*time.fs);
    if notePitch(i) ~= 0
        currFreq = midinum2frequency(notePitch(i));
    end
    freq(onset:onset+offset) = currFreq;
end
end
