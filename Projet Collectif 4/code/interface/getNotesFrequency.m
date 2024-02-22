function [freq] = getNotesFrequency(time, part)
%GETNOTESFREQUENCY - Renvoie les fréquences des notes d’une portée au cours du temps.
 % 
 % freq : vecteur contenant les fréquences des notes de la portée au cours du temps.
 % 
 % time : structure de données contenant les paramètres généraux liés au
 %     temps
 % 
 % part : matrice contenant les informations d'une portée de la partition

currFreq = 0;
freq = zeros(size(time.t));
noteOnset = [part{:,index.NOTE_ONSET}];
noteDuration = [part{:,index.NOTE_DURATION}];
notePitch = [part{:,index.MIDIPITCH}];
L = length(noteOnset);
for i=1:L
    onset = floor(noteOnset(i)*time.fs)+1;
    offset = floor(noteDuration(i)*time.fs);
    if notePitch(i) ~= 0
        currFreq = midinum2frequency(notePitch(i));
    end
    freq(onset:onset+offset) = currFreq;
end
end
