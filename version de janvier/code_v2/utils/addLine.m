function [v] = addLine(u, t1, t2, y1, y2, fs)
%returns vector u with a linear ramp between point (t1, y1) and (t2, y2)
%   à compléter...
v = u;
i1 = floor(t1*fs);
i2 = floor(t2*fs);
v(i1+1:i2) = (1:i2-i1)/fs*(y2-y1)/(t2-t1) + y1;
end

