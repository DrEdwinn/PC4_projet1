function [y] = exBow(N)
    attackSample = 1000;
    releaseSample = 10;
    endSilenceSample = 5;
    if N >= attackSample + releaseSample
        attackTemp = hamming(2*attackSample)';
        attack = attackTemp(1:attackSample);
        releaseTemp = hamming(2*(releaseSample-endSilenceSample))';
        release = [releaseTemp((releaseSample-endSilenceSample)+1:end) zeros(1,endSilenceSample)];
        y = [attack ones(1,N - attackSample - releaseSample) release];
    else
        y = hamming(N);
    end
end