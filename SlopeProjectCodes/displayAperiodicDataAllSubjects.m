% analysisChoice - 'st', 'bl' or 'combined'

% refChoice - 'none' (show raw PSDs and power) or a protocolName. For
% example, if refChoice is 'G1' then we use G1 baseline as reference.

% badTrialRejectionFlag:
% 1: Reject badElectrodes of protocolName
% 2. Reject common badElectrodes of all protocols
% 3: Reject badElectrodes of G1

% after including pairing of subjects 20-01-24

function displayAperiodicDataAllSubjects(subjectNameLists,protocolName,analysisChoice,refChoice,badEyeCondition,badTrialVersion,badElectrodeRejectionFlag,freqRangeList,axisRangeList,cutoffList,useMedianFlag,hAllPlots,pairedDataFlag)

if ~exist('protocolName','var');          protocolName='G1';            end
if ~exist('analysisChoice','var');        analysisChoice='st';          end
if ~exist('refChoice','var');             refChoice='none';             end

if ~exist('badEyeCondition','var');       badEyeCondition='ep';         end
if ~exist('badTrialVersion','var');       badTrialVersion='v8';         end
if ~exist('badElectrodeRejectionFlag','var'); badElectrodeRejectionFlag=1;  end


if ~exist('freqRangeList','var')
    freqRangeList{1} = [8 13]; % alpha
    freqRangeList{2} = [22 34]; % SG
    freqRangeList{3} = [35 65]; % FG
end
if ~exist('axisRangeList','var')
    axisRangeList{1} = [0 100];
    axisRangeList{2} = [-2.5 2.5];
    axisRangeList{3} = [-1.5 1.5];
end
if ~exist('cutoffList','var')
    cutoffList = [5 50];
end
cutoffNumElectrodes = cutoffList(1);
cutoffNumTrials = cutoffList(2);

if ~exist('useMedianFlag','var');         useMedianFlag = 0;            end
if ~exist('hAllPlots','var');             hAllPlots = [];               end

numFreqRanges = length(freqRangeList);
freqRangeColors = copper(numFreqRanges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10;
displaySettings.tickLengthMedium = [0.025 0];
displaySettings.colorNames(1,:) = [1 0 0];
displaySettings.colorNames(2,:) = [0 1 0];
titleStr{1} = 'Meditators';
titleStr{2} = 'Controls';
strList{1} = 'M';   strList{2} = 'C';

freqLims = axisRangeList{1};
yLimsPSD = axisRangeList{2};


%to change appropriately
%paramLim = {[0.5 2],[-1.5 3]};
if strcmp(refChoice,'none')
    paramLim = {[0.5 2],[-1.5 3]};
    cLimsTopo{1} = [0 3];  %exponent
    cLimsTopo{2} = [-1 2]; %offset
else
    paramLim = {[-0.5 0],[-1.5 1]};
    cLimsTopo{1} = [-0.5 0.5];  %exponent
    cLimsTopo{2} = [-1 1]; %offset
end
%nanValue = 0; % to omit nan in topoplots

%%%%%%%%%%%%%%%%%%%%%%% details and flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameterName = [{'Slope'} {'Offset '}];
refType = 'unipolar';
freqWidth = freqRangeList{1}(1);%76;
fdrFlag = 0;
bonferroniFlag = 0; %1: bonferroni correction, 2: cluster correction
if length(subjectNameLists{1}) == length(subjectNameLists{2})
    pairedFlag = 1;
else
    pairedFlag=0;
end

%%%%%%%%%%%%%%%%%%%%%%%% Get electrode groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gridType = 'EEG';
capType = 'actiCap64_UOL';
saveFolderName = 'savedData';

[electrodeGroupList,groupNameList] = getElectrodeGroups(gridType,capType);
numElecGroups = length(electrodeGroupList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(hAllPlots)
    hPSD  = getPlotHandles(1,numElecGroups,[0.05 0.55 0.6 0.3],0.02,0.02,1);
    hPower = getPlotHandles(numFreqRanges,numElecGroups,[0.05 0.05 0.6 0.45],0.02,0.02,0);
    hTopo0 = getPlotHandles(1,2,[0.675 0.7 0.3 0.15],0.02,0.02,1);
    hTopo1 = getPlotHandles(1,3,[0.675 0.55 0.3 0.13],0.02,0.02,1);
    hTopo2 = getPlotHandles(numFreqRanges,3,[0.675 0.05 0.3 0.45],0.02,0.02,1);
else
    hPSD = hAllPlots.hPSD;
    hPower = hAllPlots.hPower;
    hTopo0 = hAllPlots.hTopo0;
    hTopo1 = hAllPlots.hTopo1;
    hTopo2 = hAllPlots.hTopo2;
end

colormap("jet");
% hTopo0
montageChanlocs = showElectrodeGroups(hTopo0(1,:),capType,electrodeGroupList,groupNameList);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Protocol Position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}];
protocolPos = find(strcmp(protocolNameList,protocolName));

if ~strcmp(refChoice,'none')
    protocolPosRef = find(strcmp(protocolNameList,refChoice));
else
    protocolPosRef = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodSubjectNameLists = getGoodSubjectNameList(subjectNameLists,freqWidth,refType,badEyeCondition,badTrialVersion,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,pairedDataFlag,saveFolderName);
subjectNameLists = goodSubjectNameLists;

powerData = cell(1,2);
powerDataRef = cell(1,2);

for i=1:2
    powerDataTMP=[];    exponentDataTMP = [];   offsetDataTMP = [];
    powerDataRefTMP=[]; exponentDataTMPRef = [];   offsetDataTMPRef = [];


    for j=1:length(subjectNameLists{i})
        subjectName = subjectNameLists{i}{j};

        fileName = fullfile(pwd,saveFolderName,'FOOOF',[subjectName '_freqWidth_' num2str(freqWidth) '_' refType '__' badEyeCondition '_' badTrialVersion 'withoutKnee.mat']);
        if ~isfile(fileName)
            disp(['fileName for ' subjectName 'does not exist']);
        else
            tmpData = load(fileName);
            freqVals = tmpData.freqVals;

            [tmpPower,tmpExponent,tmpOffset] = getPowerData(tmpData,protocolPos,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials);

            if ~isempty(protocolPosRef)
                [tmpPowerRef,tmpExponentRef,tmpOffsetRef] = getPowerData(tmpData,protocolPosRef,'bl',badElectrodeRejectionFlag,cutoffNumTrials);
            end

            if isempty(protocolPosRef) % No need to worry about Ref
                if isempty(tmpPower)
                    disp(['Not enough trials for subject: ' subjectName]);
                else
                    powerDataTMP = cat(3,powerDataTMP,tmpPower);
                    exponentDataTMP = cat(2,exponentDataTMP,tmpExponent);
                    offsetDataTMP = cat(2,offsetDataTMP,tmpOffset);
                end
            else
                if isempty(tmpPowerRef) || isempty(tmpPower) % If either one is empty
                    disp(['Not enough trials for protocol or ref condition of subject: ' subjectName]);
                else
                    powerDataTMP = cat(3,powerDataTMP,tmpPower);
                    exponentDataTMP = cat(2,exponentDataTMP,tmpExponent);
                    offsetDataTMP = cat(2,offsetDataTMP,tmpOffset);

                    powerDataRefTMP = cat(3,powerDataRefTMP,tmpPowerRef);
                    exponentDataTMPRef = cat(2,exponentDataTMPRef,tmpExponentRef);
                    offsetDataTMPRef = cat(2,offsetDataTMPRef,tmpOffsetRef);

                end
            end
        end
    end
    powerData{i} = powerDataTMP;
    paramData{1}{i} = exponentDataTMP;
    paramData{2}{i} = offsetDataTMP;

    powerDataRef{i} = powerDataRefTMP;
    paramDataRef{1}{i} = exponentDataTMPRef;
    paramDataRef{2}{i} = offsetDataTMPRef;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Show Topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numElectrodes = size(powerData{1},1);
percentData = zeros(2,numElectrodes);

for i=1:2
    if isempty(protocolPosRef)
        x=powerData{i};
    else
        x=powerData{i} ./ powerDataRef{i};
    end
    numSubjects = size(x,3);

    %%%%%%%%%%%% Show percent of bad subjects per electrode %%%%%%%%%%%%%%%
    numBadSubjects = zeros(1,numElectrodes);
    for j=1:numElectrodes
        numBadSubjects(j) = sum(isnan(squeeze(x(j,1,:))));
    end
    axes(hTopo1(i)); %#ok<*LAXES>

    % Modification in the topoplot code which allows us to not interpolate across electrodes.
    percentData(i,:) = 100*(numBadSubjects/numSubjects);
    topoplot(percentData(i,:),montageChanlocs,'maplimits',[0 100],'electrodes','on','style','map','drawaxis','off','nosedir','+X'); colorbar;
    title(titleStr{i},'color',displaySettings.colorNames(i,:));
    if i==1
        ylabel('Bad subjects (%)');
    end

    %%%%%%%%%%%%%%%%%%%%%%% Show topoplots of power %%%%%%%%%%%%%%%%%%%%%%%
    for j=1:length(paramData)
        axes(hTopo2(j,i));
        if isempty(protocolPosRef)
            xParam{j}{i} = (paramData{j}{i})';
        else
            xParam{j}{i} = (paramDataRef{j}{i} - paramData{j}{i})';
        end
        if useMedianFlag
            dataTopo{j}{i} = squeeze(median(xParam{j}{i},1,'omitnan'));
        else
            dataTopo{j}{i} = squeeze(mean(xParam{j}{i},1,'omitnan'));
        end

        %dataTopo{j}{i}(find(isnan(dataTopo{j}{i}))) = nanValue;
        topoplot(dataTopo{j}{i},montageChanlocs,'electrodes','on','style','map','drawaxis','off','nosedir','+X');
        caxis(cLimsTopo{j});
        colorbar;
    end
end

%%%%%%%%%%%%%%%%%%%%%% Plot the difference of topoplots %%%%%%%%%%%%%%%%%%%
axes(hTopo1(3));
topoplot(-diff(percentData),montageChanlocs,'maplimits',[-25 25],'electrodes','on','style','map','drawaxis','off','nosedir','+X'); colorbar;

for j=1:length(paramData)
    axes(hTopo2(j,3));
    dataTopoDiff = dataTopo{j}{1} - dataTopo{j}{2};
    if ~isempty(protocolPosRef)
        dataTopoDiff = 10*dataTopoDiff;
    end
    topoplot_murty(dataTopoDiff,montageChanlocs,'electrodes','off','style','map','drawaxis','off','nosedir','+X','emarkercolors', dataTopoDiff); hold on;
    caxis([-0.5 0.5])
    %caxis([-0.2 0.2])
    sigElectrodes = findSignificantElectrodes(xParam{j},useMedianFlag,fdrFlag);
    if length(find(isnan(sigElectrodes)))~=64
        topoplot(sigElectrodes,montageChanlocs,'electrodes','on','style','blank','drawaxis','off','nosedir','+X','emarker',{'.','k',12,1});
    end
    colorbar;
end

%%%%%%%%%%%%%%%%%%%%%% Plots PSDs and power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numElecGroups
    meanPSDData = cell(1,2);
    meanPSDDataRef = cell(1,2);
    logPSDData = cell(1,2);

    badSubjectPosList = cell(1,2);
    clear pData pDataRef
    for j=1:2  % subjects
        pData{j} = powerData{j}(electrodeGroupList{i},:,:);
        numGoodElecs = length(electrodeGroupList{i}) - sum(isnan(squeeze(pData{j}(:,1,:))),1);
        badSubjectPos = find(numGoodElecs<=cutoffNumElectrodes);%union(badSubjectPos, find(numGoodElecs<=cutoffNumElectrodes));

        if ~isempty(protocolPosRef)
            pDataRef{j} = powerDataRef{j}(electrodeGroupList{i},:,:);
            numGoodElecsRef = length(electrodeGroupList{i}) - sum(isnan(squeeze(pDataRef{j}(:,1,:))),1);
            badSubjectPosRef =find(numGoodElecsRef<=cutoffNumElectrodes);% union(badSubjectPosRef,find(numGoodElecsRef<=cutoffNumElectrodes));
            badSubjectPos = union(badSubjectPos,badSubjectPosRef);%unique(cat(1,badSubjectPos,badSubjectPosRef));
        end
    end

    if pairedDataFlag
        badSubjectPosList{1} = union(badSubjectPosList{1},badSubjectPosList{2});
        badSubjectPosList{2} = badSubjectPosList{1};
    end

    for j=1:2  % subjects
        pData = powerData{j}(electrodeGroupList{i},:,:);
        for iParam=1:length(paramData);      paramDataTP{iParam}{j} = paramData{iParam}{j}(electrodeGroupList{i},:);    end

        if ~isempty(protocolPosRef)
            pDataRef = powerDataRef{j}(electrodeGroupList{i},:,:);
            for iParam=1:length(paramData);      paramDataTPRef{iParam}{j} = paramDataRef{iParam}{j}(electrodeGroupList{i},:);    end
        end
        badSubjectPos = badSubjectPosList{j};

        if ~isempty(badSubjectPos)
            disp([groupNameList{i} ', ' titleStr{j} ', '  'Not enough good electrodes for ' num2str(length(badSubjectPos)) ' subjects.']);
            pData(:,:,badSubjectPos)=[];
            for iParam=1:length(paramData)      paramDataTP{iParam}{j}(:,badSubjectPos) = [];     end
            if ~isempty(protocolPosRef)
                pDataRef(:,:,badSubjectPos)=[];
                for iParam=1:length(paramData)      paramDataTPRef{iParam}{j}(:,badSubjectPos) = [];     end
            end
        end
        meanPSDData{j} = squeeze(mean(pData,1,'omitnan'))';

        if isempty(protocolPosRef)
            logPSDData{j} = log10(meanPSDData{j});
        else
            meanPSDDataRef{j} = squeeze(mean(pDataRef,1,'omitnan'))';
            logPSDData{j} = 10*(log10(meanPSDData{j}) - log10(meanPSDDataRef{j}));
        end


        for iParam=1:length(paramData)
            if isempty(protocolPosRef)
                paramDataBar{iParam}{j} = paramDataTP{iParam}{j};
            else
                paramDataBar{iParam}{j} = paramDataTP{iParam}{j} - paramDataTPRef{iParam}{j};
            end
            meanParamData{iParam}{j} = squeeze(mean(paramDataBar{iParam}{j},1,'omitnan'));
        end

        text(30,yLimsPSD(2)-0.5*j,[titleStr{j} '(' num2str(size(meanPSDData{j},1)) ')'],'color',displaySettings.colorNames(j,:),'parent',hPSD(i));
    end
    displayAndcompareData(hPSD(i),logPSDData,freqVals,displaySettings,yLimsPSD,1,useMedianFlag,[],~pairedFlag,1,bonferroniFlag);

    title(groupNameList{i});
    xlim(hPSD(i),freqLims);

    % Add lines in PSD plots
    for jFreq=1:length(tmpData.freq_range)
        for k=1:2
            line([tmpData.freq_range{jFreq}(k) tmpData.freq_range{jFreq}(k)],yLimsPSD,'color',freqRangeColors(jFreq,:),'parent',hPSD(i));
        end
    end

    if ~isempty(protocolPosRef)
        line([0 freqVals(end)],[0 0],'color','k','parent',hPSD(i));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % bar  plots
    for iParam  =1:length(paramData)
        axes(hPower(iParam,i));
        x = categorical(strList);
        x = reordercats(x,strList);
        y=[];   err = [];
        xAll = [];   yAll = [];

        for iSub=1:2
            if useMedianFlag
                meanVal = nanmedian(meanParamData{iParam}{iSub});
                stdElectrodeGroupDataSing{iParam}{iSub}{i} = nanstd(bootstrp(10000,@nanmedian,meanParamData{iParam}{iSub}));
            else
                meanVal = nanmean(meanParamData{iParam}{iSub});
                stdElectrodeGroupDataSing{iParam}{iSub}{i} = nanstd(meanParamData{iParam}{iSub})/sqrt(length(meanParamData{iParam}{iSub}));
            end
            y = [y meanVal];
            err = [err stdElectrodeGroupDataSing{iParam}{iSub}{i}];
            yAll = [yAll meanParamData{iParam}{iSub}];
            x1 = repmat(strList{iSub},length(meanParamData{iParam}{iSub}),1)';
            if ~pairedFlag     xAll = [xAll, x1];  end
        end

        %significance level testing
        if useMedianFlag
            pD(iParam,i)=ranksum(meanParamData{iParam}{1},meanParamData{iParam}{2});
        else
            %t-test
            if pairedFlag
                [~,pD(iParam,i)]=ttest(meanParamData{iParam}{1},meanParamData{iParam}{2}); %Need to change it when the subject data stabilises
            else
                [~,pD(iParam,i)]=ttest2(meanParamData{iParam}{1},meanParamData{iParam}{2});
            end
        end


        hold on
        if pairedFlag
            for iSub = 1:2
                plot(iSub*(meanParamData{iParam}{iSub}).^0,meanParamData{iParam}{iSub},'o','MarkerFaceColor',displaySettings.colorNames(iSub,:),'MarkerEdgeColor','k');
            end

            for u= 1:min(length(meanParamData{iParam}{1}),length(meanParamData{iParam}{2}))
                plot([1 2],[meanParamData{iParam}{1}(u) meanParamData{iParam}{2}(u)],'Color', [0.5 0.5 0.5],'LineWidth',0.8);
            end
            er = errorbar([1 2],y,err,err,'LineWidth',1.5);
            xlim([0 3])
        else

            xs = categorical(cellstr(xAll'),strList);
            swarmchart(xs,yAll',20,[1 0.5 0],'filled','MarkerFaceAlpha',0.6);

            er = errorbar(x,y,err,err,'LineWidth',1.5);

        end
        er.Color =[0 0 0];
        er.LineStyle = 'none';
        er.Marker = 'diamond';
        er.MarkerSize = 10;

        pPosition = paramLim{iParam}(2)+1;

        if pD(iParam,i)<0.001
            text(1,pPosition,['p = ' num2str(pD(iParam,i),'%.1e')]);
        else
            text(1,pPosition,['p = ' num2str(round(pD(iParam,i),3),'%.3f')]);
        end

        set(hPower(iParam,i),'XTick',[1 2],'XTickLabel',[strList{1} ;strList{2}],'TickLength',[0.02, 0.02],'YTick',0:1:2,'TickDir','out','LineWidth',1,'FontWeight','bold','Box','off');
        if i==1;  ylabel(parameterName{iParam});  end

        ylim([paramLim{iParam}(1)-0.5  paramLim{iParam}(2)+1]);
        hold off
    end


end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function goodSubjectNameLists = getGoodSubjectNameList(subjectNameLists,freqWidth,refType, badEyeCondition,badTrialVersion,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,pairedDataFlag,saveFolderName)

% For unpaired case, subjects can be rejected if either data in analysis or
% ref period is bad. For paired case, a pair is rejected even if one of the two
% subjects in the pair is bad. Based on the condition, we get a new good
% subject list.

badSubjectIndex = cell(1,2);
badSubjectIndexRef = cell(1,2);

for i=1:2
    numSubjects = length(subjectNameLists{i});
    badSubjectIndexTMP = zeros(1,numSubjects);
    badSubjectIndexRefTMP = zeros(1,numSubjects);

    for j=1:numSubjects
        subjectName = subjectNameLists{i}{j};

        fileName = fullfile(pwd,saveFolderName,'FOOOF',[subjectName '_freqWidth_' num2str(freqWidth) '_' refType '__' badEyeCondition '_' badTrialVersion 'withoutKnee.mat']);

        if ~isfile(fileName)
            disp(['fileName for ' subjectName 'does not exist']);
        else
            tmpData = load(fileName);

            if isempty(getPowerData(tmpData,protocolPos,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials))
                disp(['Not enough trials for subject: ' subjectName]);
                badSubjectIndexTMP(j)=1;
            end

            if ~isempty(protocolPosRef)
                if isempty(getPowerData(tmpData,protocolPosRef,'bl',badElectrodeRejectionFlag,cutoffNumTrials))
                    disp(['Not enough trials in ref period for subject: ' subjectName]);
                    badSubjectIndexRefTMP(j)=1;
                end
            end
        end
    end

    badSubjectIndex{i} = badSubjectIndexTMP;
    badSubjectIndexRef{i} = badSubjectIndexRefTMP;
end

%%%%%%%%%%%%%%%%%%%%%%%%% Now find good subjects %%%%%%%%%%%%%%%%%%%%%%%%%%
goodSubjectNameLists = cell(1,2);

if ~pairedDataFlag
    for i=1:2
        subjectNameListTMP = subjectNameLists{i};
        if isempty(protocolPosRef)
            badPos = find(badSubjectIndex{i});
        else
            badPos = union(find(badSubjectIndex{i}),find(badSubjectIndexRef{i}));
        end
        subjectNameListTMP(badPos)=[];
        goodSubjectNameLists{i} = subjectNameListTMP;
    end
else

    if isempty(protocolPosRef)
        badPos = find(sum(cell2mat(badSubjectIndex')));
    else
        badPos = union(find(sum(cell2mat(badSubjectIndex'))),find(sum(cell2mat(badSubjectIndexRef'))));
    end

    for i=1:2
        subjectNameListTMP = subjectNameLists{i};
        subjectNameListTMP(badPos)=[];
        goodSubjectNameLists{i} = subjectNameListTMP;
    end
end
end

function [tmpPower,tmpExponent,tmpOffset] = getPowerData(tmpData,protocolPos,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials)

numTrials = tmpData.numTrials(protocolPos);
badElectrodes = getBadElectrodes(tmpData.badElectrodes,badElectrodeRejectionFlag,protocolPos);

exponentCutoff = 0.1;

if numTrials < cutoffNumTrials
    tmpPower    = [];
    tmpExponent = [];
    tmpOffset   = [];
else
    if strcmpi(analysisChoice,'st')
        tmpPower    = tmpData.SpecPower{2}{protocolPos};
        tmpExponent = tmpData.exponent{2}{protocolPos};
        tmpOffset   = tmpData.offset{2}{protocolPos};
    elseif strcmpi(analysisChoice,'bl')
        tmpPower    = tmpData.SpecPower{1}{protocolPos};
        tmpExponent = tmpData.exponent{1}{protocolPos};
        tmpOffset   = tmpData.offset{1}{protocolPos};
    else
        tmpPower = (tmpData.SpecPower{2}{protocolPos}+tmpData.SpecPower{1}{protocolPos})/2; % average
        tmpExponent = (tmpData.exponent{2}{protocolPos}+tmpData.exponent{1}{protocolPos})/2;
        tmpOffset = (tmpData.offset{2}{protocolPos}+tmpData.offset{1}{protocolPos})/2;
    end

    badExponentElec = find(tmpExponent<exponentCutoff);

    tmpPower(badElectrodes,:) = NaN;
    tmpExponent(union(badElectrodes,badExponentElec)) = NaN;
    tmpOffset(union(badElectrodes,badExponentElec)) = NaN;
end
end
function montageChanlocs = showElectrodeGroups(hPlots,capType,electrodeGroupList,groupNameList)

%%%%%%%%%%%%%%%%%%%%%% Compare with Montage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = load([capType 'Labels.mat']); montageLabels = x.montageLabels(:,2);
x = load([capType '.mat']); montageChanlocs = x.chanlocs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Topoplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(hPlots(1));
electrodeSize = 5;
numElectrodes = length(montageLabels);

for i=1:numElectrodes
    montageChanlocs(i).labels = ' ';
end

numElectrodeGroups = length(electrodeGroupList);
electrodeGroupColorList = jet(numElectrodeGroups);

for i=1:numElectrodeGroups
    topoplot(zeros(1,numElectrodes),montageChanlocs,'maplimits',[-1 1],'electrodes','on','style','map','emarker2',{electrodeGroupList{i},'o',electrodeGroupColorList(i,:),electrodeSize});
end
topoplot([],montageChanlocs,'electrodes','labels','style','blank');

axes(hPlots(2))
set(hPlots(2),'visible','off');
for i=1:numElectrodeGroups
    text(0.05,0.9-0.15*(i-1),[groupNameList{i} '(N=' num2str(length(electrodeGroupList{i})) ')'],'color',electrodeGroupColorList(i,:),'unit','normalized');
end
end
function badElectrodes = getBadElectrodes(badElectrodeList,badElectrodeRejectionFlag,protocolPos)

if badElectrodeRejectionFlag==1 % Bad electrodes for the protocol
    badElectrodes = badElectrodeList{protocolPos};
elseif badElectrodeRejectionFlag==2 % common bad electrodes for all protocols
    badElectrodes=[];
    for i=1:length(badElectrodeList)
        badElectrodes=cat(1,badElectrodes,badElectrodeList{i});
    end
    badElectrodes = unique(badElectrodes);
elseif badElectrodeRejectionFlag==3 % Bad electrodes of G1
    badElectrodes = badElectrodeList{3};
end
end


function [sigElecNum] = findSignificantElectrodes(data,useMedianFlag,fdrFlag)

if ~exist('fdrFlag','var');      fdrFlag=0;      end

allData = [];   allIDs = [];
for i=1:size(data,2)
    allData = cat(1,allData,data{i});
    allIDs = cat(1,allIDs,i+zeros(size(data{i},1),1));
end
clear p sigElecNum
sigElecNum = nan(1,size(allData,2));
for j = 1:size(allData,2)
    if useMedianFlag
        p(j)=kruskalwallis(allData(:,j),allIDs,'off');
    else
        [~,p(j)]=ttest2(data{1}(:,j),data{2}(:,j)); % only tests 2 groups
    end
end
if fdrFlag
    pAsc = sort(p,'ascend');

    for k = 1:length(pAsc)
        if ~isnan(p(k))
            p(k) = p(k)*length(p)/(find(pAsc==p(k),1));
        end
    end
end
sigElecNum(p<0.05) = 1;
end
