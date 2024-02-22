function [fPluck, fBow] = getExFunction(time, Part, ex)
%GETEXFUNCTION - Calcule les fonctions d’excitation de l’archet et/ou de pincement de la corde.
 % 
 % fPluck : vecteur contenant la fonction d’excitation pour un mode de jeu    pincé
 % 
 % fBow : vecteur contenant la fonction d’excitation pour un mode de jeu frotté
 % 
 % time : structure de données contenant les paramètres généraux liés au
 %     temps
 % 
 % part : matrice contenant les informations d'une portée de la partition
 % 
 % ex : structure de données contenant les valeurs des paramètres
 %     physiques de l'excitation

noteOnset = [Part{:,index.NOTE_ONSET}];
noteDuration = [Part{:,index.NOTE_DURATION}];
isRest = [Part{:,index.MIDIPITCH}]==0;
L = length(noteOnset);
fPluck = zeros(size(time.t));
fBow = zeros(size(time.t));

for i=1:L
    if ~isRest(i)
        onset = floor(noteOnset(i)*time.fs)+1;
        switch(ex.type)
            case 'pluck'
                offset = floor(ex.pluckDuration*time.fs);
                fPluck(onset:onset+offset-1) = hanning(offset);

            case 'bow'
                offset = floor(noteDuration(i)*time.fs);
                fBow(onset:onset+offset) = exBow(time, ex, offset+1);
                
            case 'mouette'
                offset = floor(noteDuration(i)*0.5*time.fs);
                fBow(onset:onset+offset) = exBow(time, ex, offset+1);

            case 'bowAndPluck'
                nMuffle = floor(0.1*time.fs);
                win = hanning(2*nMuffle);
                fBow(onset:onset+nMuffle-1) = win(nMuffle+1:end);
                offset = floor(ex.pluckDuration*time.fs);
                fPluck(onset:onset+offset-1) = hanning(offset);

            case 'comb'
                offset = floor(noteDuration(i)*time.fs);
                t = (1:offset)/time.fs;
                f = 20;
                a = square(2*pi*f*(2*t-t.^2/noteDuration(i)));
                fPluck(onset:onset+offset-1) = [double(a(1:end-1)>a(2:end)) 0];
        end
    end
end
end

