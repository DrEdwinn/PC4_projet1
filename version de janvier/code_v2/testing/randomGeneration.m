
Ngen = 10;

for n = 1:Ngen
    time.fs = 48000;
    parameters = randParameters();

    out = singleNote(time,parameters);

    if ~any(isnan(out))
        soundsc(out ,time.fs)
        if any(out)
            out = out/max(abs(out));
        else
            warning("PAS DE SON");
        end
    else
        warning('INPUT CONTAINS NAN VALUE WHICH CAN RESULT IN UNEXPECTED BEHAVIOR')
    end
end















