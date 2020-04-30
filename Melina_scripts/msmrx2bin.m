function iOk = msmrx2bin(dataDir, outBaseName)
iOk = -1;
impStr = 'Rhd';

smrxFiles = dir(fullfile(dataDir,'*.smrx'));
if isempty(smrxFiles)
    smrxFiles = dir(fullfile(dataDir,'*\*.smrx'));
elseif isempty(smrxFiles)
    fprintf(1,'The given directory contains no .smrx files. Please ')
    fprintf(1,'try again with another folder which do contain .smrx files.\n')
    return
end

cellSmrxFiles = {smrxFiles.name};
Nf = size(smrxFiles,1);

% Selecting files
[incFiles, iok] = listdlg('ListString',cellSmrxFiles(1,:),...
    'SelectionMode','multiple',...
    'PromptString','Select the files to join:',...
    'InitialValue',1:Nf);
if iok
    smrxFiles = smrxFiles(incFiles);
    cellSmrxFiles = {smrxFiles.name};
    Nf = size(smrxFiles,1);
else
    fprintf(1,'Cancelling...\n')
    return
end

% Merging order
fileOrder = (1:Nf)';
defInput = num2cell(num2str(fileOrder));
answr = inputdlg(cellSmrxFiles,'File order',[1, 60],defInput);
nFileOrder = str2double(answr);
nSmrxFiles = smrxFiles;
if ~isempty(answr) && sum(abs(fileOrder - nFileOrder)) ~= 0
    fprintf(1,'Changing file order...\n')
    nSmrxFiles(nFileOrder) = smrxFiles(fileOrder);
    smrxFiles = nSmrxFiles;
else
    fprintf('File order not altered\n')
end
clearvars nSmrxFiles nFileOrder
% Creating a .bin file
if ~exist('outBaseName','var') || isempty(outBaseName) ||...
        ~ischar(outBaseName)
    fprintf(1,'No outname given. Computing a name...\n')
    pathPieces = strsplit(dataDir,filesep);
    outBaseName = [pathPieces{end},'.bin'];
    fprintf(1,'File name: %s.bin\n',pathPieces{end})
end

outFullName = fullfile(dataDir,outBaseName);
if exist(outFullName,'file')
    ovwtAns = questdlg(...
        sprintf('Warning! File %s exists. Overwrite?',outFullName),...
        'Overwrite?','Yes','No','No');
    if strcmpi(ovwtAns,'no')
        fprintf('No file written.\n')
        return
    end
end
fID = fopen(outFullName,'w');
m = (2^32)/100;
for cf = 1:Nf
    cfName = fullfile(smrxFiles(cf).folder,smrxFiles(cf).name);
    FileInfo = SONXFileHeader(cfName);
    fhand = CEDS64Open(cfName);
    if fhand < 0
        fprintf(1,'The file might be opened in Spike2. Please close it and')        
        fprintf(1,' try again.\n')
        fclose(fID);
        return
    end
    totalTime = CEDS64TicksToSecs(fhand,FileInfo.maxFTime);
    
    if fhand > 0
        % Import the data.
        mxChans = CEDS64MaxChan(fhand);
        chanList = 1:mxChans;
        chTypes = zeros(mxChans,1);
        try
            heads = SONXChannelInfo(fhand,1,1);
            fch = 2;
            ch1 = 2;
            chTypes(1) = 1;
        catch
            heads = SONXChannelInfo(fhand,2,2);
            fch = 3;
            ch1 = 3;
            chTypes(2) = 1;
        end
        heads = repmat(heads,mxChans,1);
        for ch = ch1:mxChans
            chTypes(ch) = CEDS64ChanType(fhand,ch);
            if chTypes(ch) == 1
                fch = fch + 1;
                heads(ch) = SONXChannelInfo(fhand,ch,fch);
            end
        end
        heads = heads(chTypes == 1);
        chanList = chanList(chTypes == 1);
        chead = numel(heads);
        while chead >= 1
            if ~xor(isnan(str2double(heads(chead).title)),...
                    ~contains(heads(chead).title,impStr))
                heads(chead) = [];
                chanList(chead) = [];
            end
            chead = chead - 1;
        end
        multiplexerFactor = heads(1).ChanDiv;
        fs = 1 / (FileInfo.usPerTime * multiplexerFactor);
        FileInfo.SamplingFrequency = fs;
        [~, baseName, ~] = fileparts(smrxFiles(cf).name);
        save(fullfile(dataDir,...
            [baseName, '_sampling_frequency.mat']),'fs')
        display(FileInfo)
        % Determining the necessary array size to occupy approximately the
        % 75% of the available memory given that the array is int16
        memStruct = memory;
        BuffSize = 3 * memStruct.MemAvailableAllArrays / 8;
        dataPointsExp = (BuffSize / (numel(chanList) * 2));
        if heads(1).npoints < dataPointsExp
            dataPointsExp = heads(1).npoints;
        end
        wwidth = double(dataPointsExp)/fs;
        % dataPointsExp = ceil(log10(fs)+2);
        % wwidth = 10^ceil(log10(fs)+2)/fs;
        is = 1/fs;
        cw = 0;
        if abs(totalTime-(heads(1).stop - heads(1).start))
            oldTotalTime = totalTime;
            totalTime = heads(1).stop - heads(1).start;
            fprintf(1,'The length of the signals in the file seem to ')
            fprintf(1,'differ (%.3f s \\delta (%.3f ms)).\nConsidering %.3f seconds\n',...
                oldTotalTime, 1e3*(oldTotalTime - totalTime), totalTime)
        end
        while cw < totalTime
            %         if exist('Npts','var')
            %             fID = fopen(outfilename,'a');
            %         end
            % dataBuff = zeros(numel(chanList),10^dataPointsExp,'int16');
            dataBuff = zeros(numel(chanList),dataPointsExp,'int16');
            if cw <= totalTime - wwidth
                timeSegment = [cw, cw + wwidth];
            else
                timeSegment = [cw, totalTime];
                cw = totalTime * 2;
                dataBuff = zeros(numel(chanList),int32(diff(timeSegment)*fs),...
                    'int16');
            end
            shortFlag = false;
            fprintf(1,'Reading... ')
            for ch = 1:numel(chanList)
                [Npts, chanAux, ~] =...
                    SONXGetWaveformChannelSegment(fhand, chanList(ch), timeSegment,...
                    heads(ch)); %#ok<ASGLU>
                dat = int16(chanAux * m);
                try
                    dataBuff(ch,:) = dat;
                catch
                    dataBuff(ch,1:length(dat)) = dat;
                    shortFlag = true;
                end
            end
            
            fprintf(1,'done!\n')
            cw = cw + wwidth + is;
            fprintf(1,'Writting... ')
            if shortFlag
                dataBuff(:,length(dat)+1:dataPointsExp) = [];
            end
            fwrite(fID,dataBuff,'int16');
            fprintf(1,' done!\n')
            %         ftell(fID)
            %         fclose(fID);
        end
    end
    CEDS64Close(fhand);
end
fclose(fID);
fprintf('Successfully imported files!\n')
fprintf('The files merged are the following:\n')
for cf = 1:numel(smrxFiles)
    fprintf('%s\n',smrxFiles(cf).name)
end
iOk = CEDS64Close(fhand);
end