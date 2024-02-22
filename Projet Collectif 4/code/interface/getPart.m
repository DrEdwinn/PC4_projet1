function [part] = getPart(notes,partNum)
%PARTINDEX - Renvoie les notes de la portée indiquée
%
% part : un tableau cell contenant les données d’une portée, sous ensemble de score.notes
% 
% notes : un tableau cell contenant les données de toute la partition
% 
% partNum : un entier correspondant au numéro de la portée souhaitée

ind = logical([notes{:,index.PART_NUM}]'==partNum*ones(1,width(notes)));
part = notes(ind);
part = reshape(part,length(part)/width(notes),width(notes));
end
