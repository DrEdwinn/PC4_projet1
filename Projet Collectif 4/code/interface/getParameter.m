function [param] = getParameter(time, part, indParam, defaultValue)
%GETPARAMETER - Renvoie la valeur d’un paramètre au cours du temps qui peut être modifié sur un texte de portée.
 % 
 % param : vecteur contenant la valeur d’un paramètre au cours du temps
 % 
 % time : structure de données contenant les paramètres généraux liés au
 %     temps
 % 
 % part : matrice contenant les informations d'une portée de la partition
 % 
 % indParam : entier indiquant l’indice colonne du paramètre dans la  matrice score.notes
 % 
 % defaultValue : valeur initiale du paramètre au début du morceau

    param = zeros(size(time.t));
    noteOnset = [part{:,index.NOTE_ONSET}];
    targetValues = [part{:,indParam}];
    isGradual = [part{:,indParam+1}];
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
            