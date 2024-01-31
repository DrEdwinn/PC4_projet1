function [param] = getParameter(time, part, indParam, defaultValue)
    param = zeros(size(time.t));
    noteOnset = part(:,3);
    targetValues = part(:,indParam);
    isGradual = part(:,indParam+1);
    L = length(noteOnset);
    t1 = 0;
    y1 = defaultValue;
    t2 = 0;
    y2 = defaultValue;
    for n = 1:L
        if ~isnan(targetValues(n))
            t2 = noteOnset(n);
            y2 = targetValues(n);
            if isGradual(n)
                param = addLine(param, t1, t2, y1, y2, time.fs);
            else
                param = addLine(param, t1, t2, y1, y1, time.fs);
            end
            t1 = t2;
            y1 = y2;
        end
    end
    param = addLine(param, t2, time.t(end), y2, y2, time.fs);
end
            