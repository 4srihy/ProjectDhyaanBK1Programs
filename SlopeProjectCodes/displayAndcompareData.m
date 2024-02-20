function displayAndcompareData(hPlot,data,xs,displaySettings,yLims,displaySignificanceFlag,useMedianFlag,smoothSigma,nonMatchedFlag,removeDataPointFlag,bonferroniFlag,differenceFlag)

if ~exist('displaySignificanceFlag','var');       displaySignificanceFlag=0;  end
if ~exist('useMedianFlag','var');                useMedianFlag=1;             end
if ~exist('smoothSigma','var');                  smoothSigma=[];              end
if ~exist('removeDataPointFlag','var');           removeDataPointFlag=1;      end
if ~exist('bonferroniFlag','var');                bonferroniFlag = 0;         end %1 for bonferrroni correction ; 2: cluster correction
if ~exist("differenceFlag",'var') ;               differenceFlag = 0;         end

if ~exist('nonMatchedFlag','var') || isempty(nonMatchedFlag)
    nonMatchedFlag=1;
end

if useMedianFlag
    getLoc = @(g)(squeeze(median(g,1,'omitnan')));
else
    getLoc = @(g)(squeeze(mean(g,1,'omitnan')));
end

numGroups = length(data);

if ~isempty(smoothSigma)
    windowLen = 5*smoothSigma;
    window = exp(-0.5*(((1:windowLen) - (windowLen+1)/2)/smoothSigma).^2);
    window = window/sum(window); %sqrt(sum(window.^2));

    for i=1:numGroups
        data{i} = convn(data{i},window,'same');
    end
end

axes(hPlot);

if removeDataPointFlag
    removeZeroPoint=0;
    if xs(1)==0
        removeZeroPoint=1;
        xs(1)=[];
    end
    xs(length(xs))=[];
end
clear mData
for i=1:numGroups
    clear bootStat sData

    if removeDataPointFlag
        if removeZeroPoint
            data{i}(:,1)=[];
        end
        data{i}(:,size(data{i},2)) = []; %removing end point
    end
    mData{i} = getLoc(data{i});

    if ~differenceFlag

        if useMedianFlag
            bootStat = bootstrp(1000,getLoc,data{i});
            sData = std(bootStat,'omitnan');
        else
            sData = nanstd(data{i},[],1,'omitnan')/sqrt(size(data{i},1));
        end
        patch([xs';flipud(xs')],[mData{i}'-sData';flipud(mData{i}'+sData')],displaySettings.colorNames(i,:),'linestyle','none','FaceAlpha',0.4);
        hold on;
        if exist('displaySettings.lineStyle')
            plot(xs,mData{i},'color',displaySettings.colorNames(i,:),'linewidth',1.5,'LineStyle',displaySettings.lineStyle{i});
        else
            plot(xs,mData{i},'color',displaySettings.colorNames(i,:),'linewidth',1.5);
        end
        %plot(xs,mData,'linewidth',1);
    end
end

if differenceFlag
    plot(xs,mData{2}-mData{1},'color',displaySettings.colorNames(1,:),'linewidth',1.5);
end
set(gca,'fontsize',displaySettings.fontSizeLarge);
set(gca,'TickDir','out','TickLength',displaySettings.tickLengthMedium);

if exist('yLims','var') && ~isempty(yLims)
    ylim(yLims);
else
    yLims = ylim;
end

if displaySignificanceFlag % Do significance Testing %includes bonferroni and cluster correction
    pLim = 0.01;
    clusterThresh = 0.95;
    allData = [];
    allIDs = [];
    for j=1:numGroups
        allData = cat(1,allData,data{j});
        allIDs = cat(1,allIDs,j+zeros(size(data{j},1),1));
    end

    for i=1:length(xs)
        if useMedianFlag
            p(i)=kruskalwallis(allData(:,i),allIDs,'off');
        else
            if nonMatchedFlag
                [~,p(i)]=ttest2(data{1}(:,i),data{2}(:,i)); % only tests 2 groups
            else
                [~,p(i)]=ttest(data{1}(:,i),data{2}(:,i)); % only tests 2 groups
            end
        end
    end

    if bonferroniFlag == 2 %cluster based correction
        p05 = clusterCorrection(p,xs,0.05,clusterThresh);
        p01 = clusterCorrection(p,xs,0.01,clusterThresh);
    end
    % Get patch coordinates
    yVals = yLims(1)+[0 0 diff(yLims)/20 diff(yLims)/20];

    for i=1:length(xs)
        clear xMidPos xBegPos xEndPos
        xMidPos = xs(i);
        if i==1
            xBegPos = xMidPos;
        else
            xBegPos = xMidPos-(xs(i)-xs(i-1))/2;
        end
        if i==length(xs)
            xEndPos = xMidPos;
        else
            xEndPos = xMidPos+(xs(i+1)-xs(i))/2;
        end
        clear xVals; xVals = [xBegPos xEndPos xEndPos xBegPos]';

        if bonferroniFlag == 1 %Bonferroni correction
            if (p(i)<0.05/length(xs))
                patch(xVals,yVals,'g','linestyle','none');
            end
        elseif bonferroniFlag==2 %cluster correction
            if (p05(i)~=0)
                patch(xVals,yVals,'k','linestyle','none');
            end
            if (p01(i)~=0)
                patch(xVals,yVals,'g','linestyle','none');
            end
        else
            if (p(i)<0.05)
                patch(xVals,yVals,'k','linestyle','none');
            end
            if (p(i)<0.01)
                patch(xVals,yVals,'g','linestyle','none');
            end
        end
    end
end
end

function pClusterCorrected = clusterCorrection(p,xs,pLim,clusterThresh)

if ~exist('clusterThresh','var');    clusterThresh=0.95;    end
minClustSize = ceil((1-clusterThresh)*length(xs));
p1 = p;
p1((p>pLim)) = 0;
pd = bwconncomp(p1);
pClusterCorrected = zeros(1,length(xs));
if numel(pd.PixelIdxList)>0
    for iPix = 1:length(pd.PixelIdxList)
        if length(pd.PixelIdxList{iPix})>minClustSize
            pClusterCorrected(pd.PixelIdxList{iPix}) = p(pd.PixelIdxList{iPix});
        end
    end
end
end