function [f] = getExFunction(time, Part, ex)
%GETEXFUNCTION Summary of this function goes here
%   Detailed explanation goes here
f = zeros(size(time.t));
noteOnset = Part(:,3);
noteDuration = Part(:,4);
isRest = Part(:,5)==0;
L = length(noteOnset);

for i=1:L
    if ~isRest(i)
        onset = floor(noteOnset(i)*time.fs)+1;
        switch(ex.type)
            case 'pluck'
                offset = floor(ex.pluckDuration*time.fs);
                f(onset:onset+offset-1) = hanning(offset);
            case 'bow'
                offset = floor(noteDuration(i)*time.fs);
                f(onset:onset+offset) = exBow(offset+1);
            % à compléter...
        end
    end
end
end

