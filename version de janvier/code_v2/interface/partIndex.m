function [ind] = partIndex(Score,partNum)
%PARTINDEX returns logical indexing array to select the part of specified
%number
%   à compléter...
ind = Score(:,1)==partNum*ones(1,width(Score));
ind = logical(ind);
end
