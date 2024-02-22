
fs = 48000; % fréquence d'échantillonage

% lecture de la partition
Score = parseMusicXml("Hello World.musicxml");

% synthèse
audio = playScore(Score, fs);

% écoute
soundsc(audio,fs);

%%
audiowrite('out/harmo1.wav',audio2(1,:),fs);
audiowrite('out/harmo2.wav',audio2(2,:),fs);
audiowrite('out/harmo3.wav',audio2(3,:),fs);
audiowrite('out/harmo4.wav',audio2(4,:),fs);
audiowrite('out/mouette.wav',audio2(5,:),fs);
audiowrite('out/vague.wav',audio2(6,:),fs);
audiowrite('out/vague2.wav',audio2(7,:),fs);
audiowrite('out/grincements.wav',audio2(8,:),fs);
audiowrite('out/lead.wav',audio2(9,:),fs);
audiowrite('out/canon.wav',audio2(10,:),fs);
