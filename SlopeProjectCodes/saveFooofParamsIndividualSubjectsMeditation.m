%To get the fooof data for different freq ranges for each electrode

%can run for more than 1 frequency range if more than 1 centre_freq is
%specified inside the code

%%to remove freq to avoid
%note that if both the frequencies of a given peak is present, then it
%can be sustained. Only if One of the frequencies are present, then no
%peaks can be generated and fitting pose a problem. hence to be removed



%%%%%%%%getting subject details%%%%%%%%%%%

function saveFooofParamsIndividualSubjectsMeditation(folderSourceString,folderOutString,subjectName,projectName,refType,str0, knee_flag,fooofPlotFlag,rmsPlotFlag,noneElec,fooofFlag, freqRange)

if ~exist('knee_flag','var');    knee_flag=1;            end
if ~exist('fooofPlotFlag','var');    fooofPlotFlag=0;    end
if ~exist('rmsPlotFlag','var');      rmsPlotFlag = 0;    end
if ~exist('noneElec','var');     noneElec = [];          end
if ~exist('fooofFlag','var');    fooofFlag = 1;          end

if ~exist("freqRange",'var')
    allFreqRangeWidths = 86;%[106];%[76]
    centreFreq = 147;
else
    allFreqRangeWidths = freqRange(2)-freqRange(1);
    centreFreq = (freqRange(1)+freqRange(2))/2;
end

powerType = 'BLST';


switch projectName
    case 'ADGammaProject'
        protocolType = 'eyes_open';
        foldername = fullfile(folderSourceString,projectName);
        load(fullfile(foldername,protocolType,[subjectName '_' refType '_stRange_250_750.mat']));

    case 'MeditationProject'
        protocolType = 'EO1';
        foldername  = folderSourceString;
        load(fullfile(foldername,[subjectName str0 '_250_1250.mat']));

end


% SpecPower = {Type}{Protocol}(Elecs X freqVals)
SpecPower{1} = psdValsBL;
SpecPower{2} = psdValsST;

freqRes = freqVals(2)-freqVals(1);

%%%%%to give multiple frequence ranges together
minFreqVal = 4; maxFreqVal = 950;
fullFreqRange = [minFreqVal,maxFreqVal];   %overall freq_range

freqsToAvoidDeclared{1} = [0 0];
awayFromCentre = 4;
for freqAvoidcentre = 50:50:950
    freqsToAvoidDeclared{freqAvoidcentre/50 +1} = [freqAvoidcentre-awayFromCentre freqAvoidcentre+awayFromCentre];
end

for iFreqRangeWidth =1:length(allFreqRangeWidths)

    freqRangeWidth = allFreqRangeWidths(iFreqRangeWidth);%12(\pm 6),20, 24,32
    
    disp(['processing data for freqRangeWidth: ' num2str(freqRangeWidth)]);
    g=[];

    if iFreqRangeWidth==1    g=noneElec;    end

    clear freq_range freqLength
    for jFreq = 1:length(centreFreq)
        freq_range{jFreq} = [centreFreq(jFreq)-freqRangeWidth/2,centreFreq(jFreq)+freqRangeWidth/2];
        if freq_range{jFreq}(1)<minFreqVal
            freq_range{jFreq}(1) = minFreqVal;
        end

        if freq_range{jFreq}(2)>maxFreqVal
            freq_range{jFreq}(2) = maxFreqVal;
        end
    end

    numFreqRanges = length(freq_range);

    settings = struct();
    if strcmp(protocolType,'eyes_open')|| strcmp(protocolType,'eyes_closed') || strcmp(protocolType,'eyes_closed_v5')
        settings.peak_width_limits = [1,8];
    else
        settings.peak_width_limits = [4,8];
    end

    settings.max_n_peaks = 5;
    settings.min_peak_height = 0.2;
    if knee_flag      settings.aperiodic_mode = 'knee';      end
    %settings.peak_threshold = 1;

    %    g{2}{7} = [51]; % give the bad elec with none type error here
    %    g{bl/st}{protocolNum} = [electrodeNum];
    g{1}{8} = [];    g{2}{8} = [];

    for iType = 1:length(SpecPower)
        for iProt = 1:length(SpecPower{1})
            numElec = size(SpecPower{iType}{iProt},1);
            fooof_results{iType}{iProt} = [];
            exponent{iType}{iProt} = zeros(numElec,numFreqRanges);
            offset{iType}{iProt}   = zeros(numElec,numFreqRanges);
            if knee_flag    knee{iType}{iProt}     = zeros(numElec,numFreqRanges);    end
            r_SQ{iType}{iProt}     = zeros(numElec,numFreqRanges);
            icheck=[];
            icount=1;

            for ind =  1:size(SpecPower{iType}{iProt},1) %[24, 26, 29, 30, 31, 57, 58, 61, 62, 63]%   %index of the cases/control

                fooof_resultSingElec=[];
                if isnan(SpecPower{iType}{iProt} (ind,1)) || ismember(ind,g{iType}{iProt}) %|| ismember(ind,badElectrodes{iProt})
                    exponent{iType}{iProt}(ind,:) = nan; offset{iType}{iProt}(ind,:)=nan; knee{iType}{iProt}(ind,:)=nan; r_SQ{iType}{iProt}(ind,:)=nan;
                else
                    freqsToAvoid = [];
                    fooof_indFull = fooof(freqVals,SpecPower{iType}{iProt}(ind,:)',fullFreqRange,settings,true);
                    peakParams = fooof_indFull.peak_params;

                    for iPeak = 1:size(peakParams,1)
                        freqsToAvoid{iPeak} = [peakParams(iPeak,1)-peakParams(iPeak,3)/2 peakParams(iPeak,1)+peakParams(iPeak,3)/2];
                    end
                    freqsToAvoid = [freqsToAvoid freqsToAvoidDeclared];

                    if fooofFlag   % computes using fooof
                        for j=1:length(freq_range)
                            freqLength = freq_range{j}(1):freqRes:freq_range{j}(2);
                            indicesToRemove = [];
                            for iFreqAvoid = 1:length(freqsToAvoid)
                                if xor(length(find(freqLength<freqsToAvoid{iFreqAvoid}(1))),length(find(freqLength>freqsToAvoid{iFreqAvoid}(2))))
                                    indicesToRemove = [indicesToRemove find(freqLength>=freqsToAvoid{iFreqAvoid}(1) & freqLength<=freqsToAvoid{iFreqAvoid}(2))];
                                end
                            end
                            freqLength(indicesToRemove)  = nan;
                            freq_rangeNew{j} = [nanmin(freqLength), nanmax(freqLength)];

                            if (freq_rangeNew{j}(2)-freq_rangeNew{j}(1))<=4*freqRes
                                freq_rangeNew{j} = freq_range{j};
                            end

                            fooof_ind = fooof(freqVals,SpecPower{iType}{iProt} (ind,:)',freq_rangeNew{j},settings,true);
                            if knee_flag
                                [exponent{iType}{iProt}(ind,j), offset{iType}{iProt}(ind,j),knee{iType}{iProt}(ind,j),~,~,~,~,r_SQ{iType}{iProt}(ind,j)] = fooof_get_params(fooof_ind,SpecPower{iType}{iProt}(ind,:)',settings);
                            else
                                [exponent{iType}{iProt}(ind,j), offset{iType}{iProt}(ind,j),~,~,~,~,~,r_SQ{iType}{iProt}(ind,j)] = fooof_get_params(fooof_ind,SpecPower{iType}{iProt}(ind,:)',settings);
                            end

                            fooof_resultSingElec = [fooof_resultSingElec,fooof_ind];
                        end

                        fooof_results{iType}{iProt}{ind} = fooof_resultSingElec;

                    else
                        fooof_results{iType}{iProt}{ind} = getSlopesPSDBaseline_v2(log10(SpecPower{iType}{iProt}(ind,:)),freqVals,centreFreq,freqRangeWidth*ones(1,length(centreFreq)),freqsToAvoid);
                        offset{iType}{iProt}(ind,:) = log10(cell2mat(fooof_results{iType}{iProt}{ind}(1,:)));
                        exponent{iType}{iProt}(ind,:) = cell2mat(fooof_results{iType}{iProt}{ind}(2,:));
                    end

                    icheck{icount} = [ind,icount];
                    icount=icount+1;

                end
            end
            fooof_resultNumVsElecNum{iType}{iProt}=icheck;
        end
    end

    switch projectName
        case 'ADGammaProject'
            saveFolder = fullfile(folderSourceString,projectName,'FOOOF'); makeDirectory(saveFolder);

        case 'MeditationProject'
            saveFolder = fullfile(folderOutString,'savedData','FOOOF'); makeDirectory(saveFolder);
    end
    if fooofFlag
        if knee_flag
            saveFileName = fullfile(saveFolder,[subjectName '_freqWidth_' num2str(freqRangeWidth) '_' refType '_' powerType 'withKnee.mat']);
            save(saveFileName,'freqVals','SpecPower','fooof_results','exponent','offset','knee','r_SQ','fooof_resultNumVsElecNum','freqRangeWidth','centreFreq','freq_range','settings','badElectrodes','numTrials');

        else
            saveFileName = fullfile(saveFolder,[subjectName '_freqWidth_' num2str(freqRangeWidth) '_' refType '_' str0 'withoutKnee.mat']);
            save(saveFileName,'freqVals','SpecPower','fooof_results','exponent','offset','r_SQ','fooof_resultNumVsElecNum','freqRangeWidth','centreFreq','freq_range','settings','badElectrodes','numTrials');

        end

    else
        saveFileName = fullfile(saveFolder,[subjectName '_freqWidth_' num2str(freqRangeWidth) '_' refType '_' str0 'withoutFooof.mat']);
        save(saveFileName,'freqVals','SpecPower','fooof_results','exponent','offset','r_SQ','fooof_resultNumVsElecNum','freqRangeWidth','centreFreq','freq_range','badElectrodes','numTrials');
    end


    if fooofPlotFlag
        iType = 1;
        iProt = 1;% protocol type: can change if require
        %close all
        switch projectName
            case 'ADGammaProject'
                figLocation = fullfile(foldername,'FOOOF','plots'); makeDirectory(figLocation);

            case 'MeditationProject'
                figLocation = fullfile(saveFolder,'plots'); makeDirectory(figLocation);
        end
        figure(1);clf;
        for ind = 1:64
            subplot(8,8,ind);
            if ~isempty(fooof_results{iType}{iProt}{ind})
                for jFreq = 1:length(freq_range)
                    fooof_plot(fooof_results{iType}{iProt}{ind}(jFreq),false);
                end
                %xlabel(num2str(r_SQ(icheck{ind}(1))))

            end
            title(ind);
        end
        sgtitle([subjectName 'freqRangeWidth: ' num2str(freqRangeWidth)])
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,fullfile(figLocation,[subjectName 'freqRangeWidth_ ' num2str(freqRangeWidth) '_1-64_' refType '_' powerType '.fig']));

        if length(fooof_results{iType}{iProt})>64
            figure(2);clf;
            for ind = 1:48
                subplot(8,6,ind);
                if ~isempty(fooof_results{ind+64})
                    for jFreq = 1:length(freq_range)
                        fooof_plot(fooof_results{iType}{iProt}{ind+64}(jFreq),false);
                    end
                    %xlabel(num2str(r_SQ(icheck{ind}(1))))

                end
                title(ind+64);
            end
            sgtitle([subjectName ' 65-112 freqRangeWidth: ' num2str(freqRangeWidth)])
            set(gcf, 'Position', get(0, 'Screensize'));
            saveas(gcf,fullfile(figLocation,[subjectName  'freqRangeWidth_' num2str(freqRangeWidth) '_64-112_' refType '_' powerType '.fig']));
        end
    end

    if rmsPlotFlag    r_SQFull{iFreqRangeWidth} = r_SQ{iType}{iProt};     end
end

if rmsPlotFlag
    figure(3);
    xbins = 0.05:.1:.95;
    clf;
    for i=1:length(centreFreq)
        subplot(4,4,i)
        r_SQPlot = [];
        for iFreqRange = 1:length(allFreqRangeWidths)
            r_SQPlot = [r_SQPlot , r_SQFull{iFreqRange}(:,i)];
        end
        hist(r_SQPlot,xbins); hold on;
        title(centreFreq(i));
        xlim([0.4 1]);
    end
    %legend(mat2str(allFreqRangeWidths));
    sgtitle( [subjectName ' R^2 freqRangeWidths: ' mat2str(allFreqRangeWidths)]);
    figLocation = fullfile(foldername,'FOOOF','rmsPlots'); makeDirectory(figLocation);
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,fullfile(figLocation,[subjectName '_' refType '_' powerType '.fig']));
end

end

