function [y] = normalise(x)
    if max(abs(x)) ~= 0
        y = x/max(abs(x));
    else 
        y = x;
    end
end