prompt = "commandes : \n 0 : quitter \n 1 : passer au son suivant \n 2 : réécouter \n 3 : enregistrer le preset\n";
in = -1;
time.fs = 48000;
out = 0;
while in ~=0
    clc
    while any(isnan(out)) || ~any(out)
        parameters = randParameters();
        out = singleNote(time,parameters);
    end
    
    soundsc(out,time.fs);
    in = input(prompt);

    switch(in)
        case 1
            out = zeros(size(out));
        case 2
            soundsc(out,time.fs);
        case 3
            name = input("entrez le nom du preset\n","s");
            savePreset(parameters,name);
    end
end

