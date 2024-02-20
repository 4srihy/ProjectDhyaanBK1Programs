clear; close all

projectName = 'MeditationProject';
[allSubjectNames,expDateList] = getDemographicDetails('BK1');
[goodSubjectList, meditatorList, controlList] = getGoodSubjectsBK1;
folderSourceString = 'C:\Users\Srishty\OneDrive - Indian Institute of Science\Documents\supratim\meditationProject\codes\ProjectDhyaanBK1Programs\powerProjectCodes\savedData';
%folderOutString  = 'C:\Users\Srishty\OneDrive - Indian Institute of Science\Documents\supratim\meditationProject\codes\ProjectDhyaanBK1Programs\SlopeProjectCodes';
folderOutString  = 'C:\Users\Srishty\OneDrive - Indian Institute of Science';

refType = 'unipolar';
freqRange = [104, 190]; %ensure that centreFreq is integer
fooofFlag = 1; %if 0, would calculate slope using fit method
knee_flag = 0;
fooofPlotFlag = 0;
rmsPlotFlag = 0;
noneElec = []; % for some individual subjects, fooof fitting gives nonetype error. you have to identify the electrode from the code and the reject it

badEyeCondition = 'ep'; % use 'wo' for without, 'ep' for eye position and 'em' for eye movement
badTrialVersion = 'v8';
badTrialNameStr = ['_' badEyeCondition '_' badTrialVersion];

useTheseIndices = 1:length(goodSubjectList);

for i=setdiff(1:35,[12,17,37,40,71])
    subjectName = goodSubjectList{i};
    disp([num2str(i) ':' subjectName]);
    saveFooofParamsIndividualSubjectsMeditation(folderSourceString,folderOutString,subjectName,projectName,refType,badTrialNameStr, knee_flag,fooofPlotFlag,rmsPlotFlag,noneElec,fooofFlag,freqRange);
end

%% None type object error in electrodes: getNan for bad elecs
%ep i=12-053DR  g{1}{7}(33) , 17-044PN: g{1}{6}(49) , 37-021PB :g{2}{7}(17) , 40-027SM:g{2}{6}(52) , 74-102AS: g{2}{7}(51)
% wo i=12 g{1}{7}(33)   ,37: g{2}{7}(17), 40: g{2}{6}(19),  61: g{1}{5}(22),    74:g{2}{7}(51)

