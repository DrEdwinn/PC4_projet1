
function [y] = positivePart(x)
    if x > 0
        y = x;
    else
        y = 0;
    end
end
