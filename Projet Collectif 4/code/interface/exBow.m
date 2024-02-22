function [y] = exBow(time, ex, exDurationSample)
    attackSample = floor(ex.bowAttack * time.fs);
    releaseSample = 500;
    endSilenceSample = 500;
    
    y = ones(1,exDurationSample);
    win1 = hanning(2*attackSample);
    win2 = hanning(2*releaseSample);
    y(1:attackSample) = win1(1:attackSample);
    releaseOnset = exDurationSample - endSilenceSample - releaseSample;
    y(releaseOnset:exDurationSample-endSilenceSample) = win2(end-releaseSample:end);
    y(exDurationSample - endSilenceSample + 1:end) = 0;
    y = y + ex.bowModAmount * (sin(2 * pi * ex.bowModFreq * (0:exDurationSample-1)/time.fs)+1)/2;
end