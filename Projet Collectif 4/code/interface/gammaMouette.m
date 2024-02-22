function [gamma] = gammaMouette(time, part, defaultGamma)
%GAMMAMOUETTE - Renvoie la valeur du paramètre gamma pour obtenir un son de mouette.
 % 
 % gamma : vecteur contenant la valeur du paramètre gamma au cours du temps
 % 
 % time : structure de données contenant les paramètres généraux liés au
 %     temps
 % 
 % part : matrice contenant les informations d'une portée de la partition
 % 
 % defaultGamma : la valeur initiale de gamma pour cette corde

    noteOnset = [part{:,index.NOTE_ONSET}];
    noteDuration = [part{:,index.NOTE_DURATION}];
    noteDuration = floor(noteDuration*100)/100;
    notePitch = [part{:,index.MIDIPITCH}];
    L = length(noteOnset);
    
    gamma = zeros(size(time.t));
    for i=1:L
        if notePitch(i)~=0
            onset = floor(noteOnset(i)*time.fs)+1;
            offset = floor(noteDuration(i)*time.fs);

            tq = [0 .3 .6 .7 .95 1]*noteDuration(i); % dessin de l'enveloppe gamma souhaitée
            gammaq = [0 .6 .4 1 .05 0];
            t = (0:.01:noteDuration(i));
  
            gammaNote = bezier.eval(interp1(tq,gammaq,t)',noteDuration(i)*time.fs); % utilisation de courbe de bézier pour arrondir l'enveloppe
            gammaNote = gammaNote/max(gammaNote);
            gamma(onset:onset+offset-1) = gammaNote*500;
        end
    end
    gamma = gamma + defaultGamma;
end

