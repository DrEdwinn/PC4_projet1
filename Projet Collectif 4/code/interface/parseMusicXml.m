%PARSEMUSICXML - Extrait les données depuis un fichier au format .musicxml
% 
% score : structure de donnée comportant 3 champs :
% score.notes : un tableau cell où est stocké toutes les donnée de la partition. Ce tableau comporte autant de lignes que le nombre total de notes présentes dans la partition et 16 colonnes dont chacune correspond à une information, comme décrit ci dessous.
% 1 : nom des presets si indiqués sur la partition
% 2 : le numéro de la portée
% 3 : le numéro de la mesure
% 4 : le début de la note, en s
% 5 : la durée de la note, en s
% 6 : le pitch midi de la note
% 7 : valeur de detune si indiqué sur la partition
% 8 : changement progressif de detune indiqué sur la partition ? (0/1)
% 9 : valeur de kappa si indiqué sur la partition
% 10 : changement progressif de kappa indiqué sur la partition ? (0/1)
% 11 : valeur de position de ressort (xr) si indiqué sur la partition
% 12 : changement progressif de xr indiqué sur la partition ? (0/1)
% 13 : valeur de raideur du ressort linéaire (omega_r0) si indiqué sur la partition
% 14 : changement progressif de omega_r0 ? (0/1)
% 15 : valeur de raideur du ressort cubique (omega_r1) si indiqué sur la partition
% 16 : changement progressif de omega_r1 ? (0/1)
% 
% score.duration : la durée du morceau en secondes
% 
% score.partCount : le nombre de portées dans la partition
% 
% filename : chaîne de caractères indiquant le nom du fichier à lire
%
% Copyright 2018 Fabio Jose Muneratti Ortega.
%
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
% DEALINGS IN THE SOFTWARE.
%
function [score] = parseMusicXml(filename)

    % Create the DocumentBuilder
    builder = javaMethod('newInstance', 'javax.xml.parsers.DocumentBuilderFactory');
    
    % Disable validation (because of MATLAB's xmlread bug)
    builder.setFeature('http://apache.org/xml/features/nonvalidating/load-external-dtd', false);
    
    dom = xmlread(filename, builder);
    
    ind = 1; % note array index 
    % parse first 'part'
    parts = dom.getElementsByTagName('part');
    measures = parts.item(0).getElementsByTagName('measure');
    mxml = cell(measures.getLength,16); % pre-allocation
    
    currKey = 0; % default C maj
    currTempo = 100; % default tempo 100 bpm
    
    partCount = parts.getLength();
    for iPart = 0:partCount-1
        part = parts.item(iPart);
        partNum = double(erase(string(part.getAttributes().getNamedItem('id').getValue()),'P'));
        measures = part.getElementsByTagName('measure');
        
        beatCounter = 0;
        currBeats = NaN;
        currOrnam = 'n'; % no ornamentation
        detuneValue = NaN ; 
        detuneGradualChange = 0;
        kappaValue = NaN ; 
        kappaGradualChange = 0;
        xrValue = NaN;
        xrGradualChange = 0;
        omr0Value = NaN;
        omr0GradualChange = 0;
        omr1Value = NaN;
        omr1GradualChange = 0;
        presetName = '';
        slurStarted = 0;
        
        for iMeasure = 0:(measures.getLength()-1)
            % parse each measure
            measure = measures.item(iMeasure);
            measureNum = javaMethod('parseInt', 'java.lang.Integer', measure.getAttributes(). ...
                getNamedItem('number').getValue());
            % reset sharps/flats for new measure;
            currAccidentals = accidentalsForKey(currKey);
            % parse each measure element
            for iNode = 0:(measure.getLength()-1)
                switch char(measure.item(iNode).getNodeName())
                    case 'attributes'
                        [currKey, currBeats] = parseAttributes(measure.item(iNode));
                        if isnan(currKey)
                            error("nan fifths")
                        end
                        currAccidentals = accidentalsForKey(currKey);
                    case 'note'
                        [st, oc, d, a, v, o, s, ac] = parseNote(measure.item(iNode));
                        if strcmp(char(d),'measureDuration')
                            d = currBeats;
                        end
                        if (~isnan(ac))
                            % Note indicated a change of accidentals.
                            currAccidentals = updateAccidentals(st, ac, currAccidentals);
                        end
                        if (o ~= 'g') % if the note isn't a grace note, include it in output
                            p = midiPitch(st, oc, currAccidentals);
                            if (o == '+' && a == 'l')
                                % if this note is linked with a tie, change
                                % duration of the previous note and make
                                % this note a rest
                                mxml{ind-1, index.NOTE_DURATION} = mxml{ind-1, index.NOTE_DURATION} + d*60.0/currTempo;
                                p = 0;
                            end
                            if currOrnam == 'n' % if there is no grace note to indicate
                                currOrnam = o; % indicated ornamentation is as parsed
                            end
                            mxml(ind, index.PRESET_NAME) = {presetName};
                            mxml(ind, index.PART_NUM) = {partNum};
                            mxml(ind, index.MEASURE_NUM) = {measureNum};
                            mxml(ind, index.NOTE_ONSET) = {beatCounter*60.0/currTempo};
                            mxml(ind, index.NOTE_DURATION) = {d*60.0/currTempo};
                            mxml(ind, index.MIDIPITCH) = {p};
                            mxml(ind, index.DETUNE_VALUE) = {detuneValue};
                            mxml(ind, index.DETUNE_ISGRADUAL) = {detuneGradualChange};
                            mxml(ind, index.KAPPA_VALUE) = {kappaValue};
                            mxml(ind, index.KAPPA_ISGRADUAL) = {kappaGradualChange};
                            mxml(ind, index.XR_VALUE) = {xrValue};
                            mxml(ind, index.XR_ISGRADUAL) = {xrGradualChange};
                            mxml(ind, index.OMR0_VALUE) = {omr0Value};
                            mxml(ind, index.OMR0_ISGRADUAL) = {omr0GradualChange};
                            mxml(ind, index.OMR1_VALUE) = {omr1Value};
                            mxml(ind, index.OMR1_ISGRADUAL) = {omr1GradualChange};
                            
                            detuneValue = NaN;
                            kappaValue = NaN;
                            xrValue = NaN;
                            omr0Value = NaN;
                            omr1Value = NaN;
                            presetName = '';
                            currOrnam = 'n';
                            ind = ind + 1;
                            
                            beatCounter = beatCounter + d;
                            slurStarted = s + slurStarted;
                        else
                            currOrnam = o;
                        end
                    case 'direction'
                        typeList = measure.item(iNode).getElementsByTagName('direction-type').item(0).getChildNodes();
                        for k = 0:(typeList.getLength()-1)
                            type = typeList.item(k);
                            switch char(type.getNodeName())
                                case 'metronome'
                                    unit = char(type.getElementsByTagName('beat-unit'). ...
                                        item(0).getFirstChild().getNodeValue());
                                    tempo = javaMethod('parseInt', 'java.lang.Integer', ...
                                        type.getElementsByTagName('per-minute').item(0).getFirstChild().getNodeValue());
                                    if (strcmp(unit, 'quarter'))
                                        currTempo = tempo;
                                    elseif (strcmp(unit, 'eighth'))
                                        currTempo = tempo/2;
                                    else
                                        currTempo = tempo;
                                    end
                                case 'dynamics'
                                    dynNodes = type.getChildNodes();
                                    currDyn = '0'; % unknown
                                    for l = 0:(dynNodes.getLength()-1)
                                        dyn = char(dynNodes.item(l).getNodeName());
                                        if (strcmpi(dyn, 'fz') || strcmpi(dyn, 'sfz') || strcmpi(dyn, 'sf'))
                                            currDyn = 's';
                                        elseif (strcmpi(dyn, 'fff'))
                                            currDyn = '9';
                                        elseif (strcmpi(dyn, 'ff'))
                                            currDyn = '8';
                                        elseif (strcmpi(dyn, 'f'))
                                            currDyn = '7';
                                        elseif (strcmpi(dyn, 'mf'))
                                            currDyn = '6';
                                        elseif (strcmpi(dyn, 'mp'))
                                            currDyn = '5';
                                        elseif (strcmpi(dyn, 'p'))
                                            currDyn = '4';
                                        elseif (strcmpi(dyn, 'pp'))
                                            currDyn = '3';
                                        elseif (strcmpi(dyn, 'ppp'))
                                            currDyn = '2';
                                        end
                                    end
                                    currDynOnset = beatCounter;
                                case 'wedge'
                                    wt = char(type.getAttributes().getNamedItem('type').getValue());
                                    switch wt
                                        case 'stop'
                                            currWedge = 'n';
                                        otherwise
                                            currWedge = wt(1);
                                    end
                                case 'words'
                                    if ~isempty(type.item(0))
                                        word = char(type.item(0).getNodeValue());
                                        [pn, tv, g] = parseWord(word);
                                        switch pn
                                            case {'detune','d','dt'}
                                                detuneValue = str2double(tv);
                                                detuneGradualChange = g;
                                                
                                            case {'kappa','k'}
                                                kappaValue = str2double(tv);
                                                kappaGradualChange = g;
                                            case 'xr'
                                                xrValue = str2double(tv);
                                                xrGradualChange = g;
                                            case {'omr0','omega_r0'}
                                                omr0Value = str2double(tv);
                                                omr0GradualChange = g;
                                            case {'omr1','omega_r1'}
                                                omr1Value = str2double(tv);
                                                omr1GradualChange = g;
                                            case 'preset'
                                                presetName = tv;
                                        end
                                    end
                            end
                        end
                end
            end
        end
    end
    score.notes = mxml;
    score.duration = mxml{end,index.NOTE_ONSET} + mxml{end,index.NOTE_DURATION};
    score.partCount = partCount;
end

function accid = accidentalsForKey(fifths)
    fifthsMat = [0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 1, 0, 0, 0;
        1, 0, 0, 1, 0, 0, 0;
        1, 0, 0, 1, 1, 0, 0;
        1, 1, 0, 1, 1, 0, 0;
        1, 1, 0, 1, 1, 1, 0;
        1, 1, 1, 1, 1, 1, 0;
        1, 1, 1, 1, 1, 1, 1];
    
    if (fifths >= 0)
        accid = fifthsMat(fifths+1, :);
    else
        accid = fifthsMat(8+fifths,:) - 1;
    end
end

function acList = updateAccidentals(step, newAc, acList)
    if strcmpi('C', step)
        acList(1) = newAc;
    elseif strcmpi('D', step)
        acList(2) = newAc;
    elseif strcmpi('E', step)
        acList(3) = newAc;
    elseif strcmpi('F', step)
        acList(4) = newAc;
    elseif strcmpi('G', step)
        acList(5) = newAc;
    elseif strcmpi('A', step)
        acList(6) = newAc;
    elseif strcmpi('B', step)
        acList(7) = newAc;
    end
end

function p = midiPitch(step, octave, accid)
p = 12 + octave*12;
if strcmpi('C', step)
    p = p + accid(1);
elseif strcmpi('D', step)
    p = p + 2 + accid(2);
elseif strcmpi('E', step)
    p = p + 4 + accid(3);
elseif strcmpi('F', step)
    p = p + 5 + accid(4);
elseif strcmpi('G', step)
    p = p + 7 + accid(5);
elseif strcmpi('A', step)
    p = p + 9 + accid(6);
elseif strcmpi('B', step)
    p = p + 11 + accid(7);
else %rest
    p = 0;
end
end

function [step, octave, duration, art, vib, orn, slur, accid] = parseNote(note)
    nodes = note.getChildNodes();
    hasDot = 0;
    art = 'l'; % default: legato
    vib = 0;
    orn = 'n'; % default: no ornamentation
    slur = 0;
    step = NaN;
    accid = NaN;
    octave = 0;
    duration = 0;
    timeModNum = 1;
    timeModDen = 1;
    for i = 0:(nodes.getLength()-1)
        elmt = nodes.item(i);
        switch char(elmt.getNodeName())
            case 'pitch'
                step = char(elmt.getElementsByTagName( ...
                'step').item(0).getFirstChild().getNodeValue());
                octave = javaMethod('parseInt', 'java.lang.Integer', elmt.getElementsByTagName( ...
                'octave').item(0).getFirstChild().getNodeValue());
            case 'dot'
                hasDot = 1;
            case 'rest'
                step = 'rest';
                if (elmt.hasAttributes()) % if it's a rest measure
                    duration = 'measureDuration';
                end
            case 'tie'
                if (strcmpi('stop', ...
                        char(elmt.getAttributes().getNamedItem('type'). ...
                        getValue())))
                    orn = '+';
                end
            case 'type'
                switch char(elmt.getFirstChild().getNodeValue())
                    case 'whole'
                        duration = 4;
                    case 'half'
                        duration = 2;
                    case 'quarter'
                        duration = 1;
                    case 'eighth'
                        duration = 0.5;
                    case '16th'
                        duration = 0.25;
                    case '32nd'
                        duration = 0.125;
                end
            case 'time-modification'
                timeModNum = javaMethod('parseInt', 'java.lang.Integer', elmt.getElementsByTagName( ...
                'normal-notes').item(0).getFirstChild().getNodeValue());
                timeModDen = javaMethod('parseInt', 'java.lang.Integer', elmt.getElementsByTagName( ...
                'actual-notes').item(0).getFirstChild().getNodeValue());
            case 'accidental'
                switch char(elmt.getFirstChild().getNodeValue())
                    case 'sharp'
                        accid = 1;
                    case 'flat'
                        accid = -1;
                    case 'double-sharp'
                        accid = 2;
                    case 'double-flat'
                        accid = -2;
                    case 'natural'
                        accid = 0;
                end
            case 'grace'
                orn = 'g'; % todo: accacciatura or appogiatura
            case 'notations'
                nots = elmt.getChildNodes();
                for j = 0:(nots.getLength()-1)
                    switch char(nots.item(j).getNodeName())
                        case 'slur'
                            t = nots.item(j).getAttributes().getNamedItem( ...
                                'type').getValue();
                            if strcmpi(t, 'start')
                                slur = 1;
                            else % stop
                                slur = -1;
                            end
                        case 'articulations'
                            art = nots.item(j).getChildNodes();
                            if (art.getElementsByTagName('accent').getLength() > 0)
                                art = '<';
                            elseif (art.getElementsByTagName('staccato').getLength() > 0)
                                art = '.';
                            elseif (art.getElementsByTagName('tenuto').getLength() > 0)
                                art = '-';
                            else
                                art = 'l';
                            end
                        case 'technical'
                            % todo
                        case 'ornaments'
                            ornmt = nots.item(j).getChildNodes();
                            if (ornmt.getElementsByTagName('inverted-mordent').getLength() > 0)
                                orn = 't';
                            end
                    end
                end    
        end   
    end
    duration = duration *timeModNum/timeModDen;
    if hasDot
        duration = duration * 1.5;
    end
end

function [key, beat, beatType] = parseAttributes(attr)

    key = NaN;
    beat = NaN;
    beatType = NaN;

    % change of key
    k = attr.getElementsByTagName('key');
    if (k.getLength() > 0)
        fifths = k.item(0).getElementsByTagName('fifths');
        if (fifths.getLength() > 0)
            key = javaMethod('parseInt', 'java.lang.Integer', fifths.item(0).getFirstChild().getNodeValue());
        end
    end
    
    % change of time signature
    time = attr.getElementsByTagName('time');
    if (time.getLength() > 0)
        beats = time.item(0).getElementsByTagName('beats');
        if (beats.getLength() > 0)
            beat = javaMethod('parseInt', 'java.lang.Integer', beats.item(0).getFirstChild().getNodeValue());
        end
        beatType = time.item(0).getElementsByTagName('beat-type');
        if (beatType.getLength() > 0)
            beatType = javaMethod('parseInt', 'java.lang.Integer', beatType.item(0).getFirstChild().getNodeValue());
        end
    end
    
    % ignoring clef...
end

function [paramName, targetValue, isGradual] = parseWord(s)
    paramName = "";
    targetValue = NaN;
    isGradual = 0;
    s = erase(s,' '); % erasing whitespaces
    isExpressionValid = 1;
    if contains(s, '=')
        out = split(s,'=');
        isGradual = 0;
    elseif contains(s,'~')
        out = split(s,'~');
        isGradual = 1;
    else
        isExpressionValid = 0;
    end
    

    if (isExpressionValid)
        paramName = out{1,1};
        targetValue = out{2,1};
    end
end

