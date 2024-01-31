% return corresponding frequency of given midi pitch
function [f] = midinum2frequency(midinum)
f = 8.1757989156*2.^((midinum)/12);
end