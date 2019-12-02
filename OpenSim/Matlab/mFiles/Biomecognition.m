%% Biomecognition Paper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hossein Mokhtarzadeh
% Last Date: Oct 7th 2019
% Initial Date: Sep 5th 2019
% What does it do?

% This paper here: ~\Projects\Publications\Biomecognition_PaperMethod so Biomcognition_DualTask_HM_v12
% I am focusng on 1 person SH and get all the results
% 1. a table for peak GRFs, major muscles, mean Mos, gait param, cost, performace, 

% 2. Focuso on these activity only no pertub and 
% 3. make average of all 4 optinos when avaibale. (7 activies and 4 times each, n=1)
trials = {'S+T01','S+T02','S+T03','S+T04', ...
'W101','W102','W103','W104',...
    'W1+T01','W1+T02','W1+T03','W1+T04',...
    'W201','W202','W203','W204',...
    'W2+T01','W2+T02','W2+T03','W2+T04',...
    'W301','W302','W303','W304',...
    'W3+T01','W3+T02','W3+T03','W3+T04',};

% This program is to develop and create results for biomecognition paper e.g. 
% GRF, Joint Kinematisc, etc using Visualizegetbasicmot.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% source folder
clear 
clc
close all
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH';
cd(source)
mydir  = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1);
tr = fullfile(newdir,'Trials', 'trials.mat'); % create trials.mat in a folder called Trials
 %these are
  trials = {'S+T01','S+T02','S+T03','S+T04', ...
'W101','W102','W103','W104','W1+R+ACC','W1+R+DEC','W1+L+ACC','W1+L+DEC',...
    'W1+T01','W1+T02','W1+T03','W1+T04','W1+T+R+ACC','W1+T+R+DEC','W1+T+L+ACC','W1+T+L+DEC',...
    'W201','W202','W203','W204','W2+R+ACC','W2+R+DEC','W2+L+ACC','W2+L+DEC',...
    'W2+T01','W2+T02','W2+T03','W2+T04','W2+T+R+ACC','W2+T+R+DEC','W2+T+L+ACC','W2+T+L+DEC', ...
    'W301','W302','W303','W304','W3+R+ACC','W3+R+DEC','W3+L+ACC','W3+L+DEC',...
    'W3+T01','W3+T02','W3+T03','W3+T04','W3+T+R+ACC','W3+T+R+DEC','W3+T+L+ACC','W3+T+L+DEC'};
%   trials = {'S+T01' (1),'S+T02'(2),'S+T03'(3),'S+T04'(4), ...
% 'W101'(5),'W102'(6),'W103'(7),'W104'(8),'W1+R+ACC'(9),'W1+R+DEC'(10),'W1+L+ACC'(11),'W1+L+DEC'(12),...
%     'W1+T01'(13),'W1+T02'(14),'W1+T03'(15),'W1+T04'(16),'W1+T+R+ACC'(17),'W1+T+R+DEC'(18),'W1+T+L+ACC'(19),'W1+T+L+DEC'(20),...
%     'W201'(21),'W202'(22),'W203'(23),'W204'(24),'W2+R+ACC'(25),'W2+R+DEC'(26),'W2+L+ACC'(27),'W2+L+DEC'(28),...
%     'W2+T01'(29),'W2+T02'(30),'W2+T03'(31),'W2+T04'(32),'W2+T+R+ACC'(33),'W2+T+R+DEC'(34),'W2+T+L+ACC'(35),'W2+T+L+DEC'(36), ...
%     'W301'(37),'W302'(38),'W303'(39),'W304'(40),'W3+R+ACC'(41),'W3+R+DEC'(42),'W3+L+ACC'(43),'W3+L+DEC'(44),...
%     'W3+T01'(45),'W3+T02'(46),'W3+T03'(47),'W3+T04'(48),'W3+T+R+ACC'(49),'W3+T+R+DEC'(50),'W3+T+L+ACC'(51),'W3+T+L+DEC'(52)};


% load(tr)
grf = fullfile(source, 'grfResults'); trialgrf = string(strcat(trials{10},'.mot'));
trc = fullfile(source, 'trcResults');trialtrc = string(strcat(trials{10},'.trc'));%'W304.trc'
SO = fullfile(source, 'SOResults'); trialSO = string(strcat(trials{19},'_StaticOptimization_force.sto'));
ID = fullfile(source, 'IDResults');trialID = string(strcat(trials{10},'_ID.sto'));
IK = fullfile(source, 'IKResults');trialIK = string(strcat(trials{10},'_IK.mot'));

%% Biomechanical variables

%% Table for max, min, stat, for the following variables in all 7 second 4 intervals,first/last cycle:
% TODos:
% 1. do the table correclty for frist/last average 4 times, etc (maybe not)
% 2. Do the table also for whole 7 seconds and 4 times (this is first) all
% 7secod is better be consisteant for graph, and table...
% 3.  make three tables (muscles, kinemtic, kineti, and then cost of them)
% 4. then think of what graphs to show. maybe graob first/last is good
% 4. no stat needed as n=1, just show how much % increase decrease in
% loading, etc cost happenes e.e.g cost 20% increase etc

% varibales for the table: {GRFs, GRFs, GRFv, timing, GRFimpact, maxKnee, maxhip,maxankle, ...
% moments (hip, knee, anke), MoS, Mos, stepL,StepFre, Cost, perforance, .. }
% S: Static, W: Walking, T: Target
% variable/Trial | S | W1 | W1+T | W2 | W2+T | W3 | W3+T | Stat |
% GRF            | S | W1 | W1+T | W2 | W2+T | W3 | W3+T | Stat |
% MoS            | S | W1 | W1+T | W2 | W2+T | W3 | W3+T | Stat |
% Muslces        | S | W1 | W1+T | W2 | W2+T | W3 | W3+T | Stat |
% Moments        | S | W1 | W1+T | W2 | W2+T | W3 | W3+T | Stat |
% Motions        | S | W1 | W1+T | W2 | W2+T | W3 | W3+T | Stat |
% GaitVar        | S | W1 | W1+T | W2 | W2+T | W3 | W3+T | Stat |
% ...            | S | W1 | W1+T | W2 | W2+T | W3 | W3+T | Stat |

% vars = {'vGRF'}; % then add any variable from the table

% 0. Path and data location etc and whic trials
clear;clc;close all
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';
load(fullfile(source, 'DFLOW\dflow.mat'))
load(fullfile(source, 'trcResults\gaitEvents.mat'))
trials = {'S+T01','S+T02','S+T03','S+T04', ...
'W101','W102','W103','W104',...
    'W1+T01','W1+T02','W1+T03','W1+T04',...
    'W201','W202','W203','W204',...
    'W2+T01','W2+T02','W2+T03','W2+T04',...
    'W301','W302','W303','W304',...
    'W3+T01','W3+T02','W3+T03','W3+T04'};

% load all the data GRF, muscles, mos, etc
load(fullfile(source,'grfResults\GRFs.mat'))
load(fullfile(source,'SOResults\SOs.mat'))
load(fullfile(source,'IDResults\IDs.mat'))
load(fullfile(source,'IKResults\IKs.mat'))
load(fullfile(source,'trcResults\gaitEvents.mat'))
load(fullfile(source,'trcResults\gaitUtils.mat'))

trialsUtils= gaitUtils.trialNames;
fields = fieldnames(GRFs(6).data);
for i = 1:length(GRFs)
    trials2{i}= GRFs(i).name;
   trialsEvent{i}= gaitEvents(i).name;
   
end

for i = 1:length(IDs)
 
    trialsID{i}= extractBefore(IDs(i).name,'_');
    
end


for i = 1:length(IKs)
   
    trialsIK{i}= extractBefore(IKs(i).name,'_');
   
end


for i = 1:length(SOs)
   
    try 
    
    trialsSO{i}=  extractBefore(SOs(i).name,'_');
    
    catch 
        fprintf('Inconsistent data in iteration %s, skipped.\n', i);
    end
end

for j=1:length(trials)
% 1. load vars{1} e.g. vGRF load GRFs

% load(fullfile(source,'grfResults\GRFs.mat'))



[ind, missTrials]=trialIndexFinder(string(trials{j}),trials2);
%get gaitEvnts trials from all
[indEvent, ~]=trialIndexFinder(string(trials{j}),trialsEvent);
%get gaitutils trials from all
[indUtils, ~]=trialIndexFinder(string(trials{j}),trialsEvent);
try 
[indSO, missTrialsSO]=trialIndexFinder(string(trials{j}),trialsSO);
catch 
    print('chcek later')
end
try
[indID, missTrialsID]=trialIndexFinder(string(trials{j}),trialsID);
catch 
    print('chcek later')
end
try
[indIK, missTrialsIK]=trialIndexFinder(string(trials{j}),trialsIK);
catch 
    print('chcek later')
end

% ini = dflow.Perturb(28).iniPerturb;
% last = dflow.Perturb(28).endPertub;


% 2. find first and last cycle timing , you may code this, done it before
%lets find first and last cycles
ini1 = gaitEvents(ind).data.HSLocL(1);
ini2 = gaitEvents(ind).data.HSLocL(2);
last1 = gaitEvents(ind).data.HSLocL(end-1);
last2 = gaitEvents(ind).data.HSLocL(end);
[~,locIni1]=min(abs(ini1-GRFs(ind).data.(fields{end})));
[~,locIni2]=min(abs(ini2-GRFs(ind).data.(fields{end})));
[~,locLast1]=min(abs(last1-GRFs(ind).data.(fields{end})));
[~,locLast2]=min(abs(last2-GRFs(ind).data.(fields{end})));
subplot(2,1,1);plot(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{2})(locIni1:locIni2));
subplot(2,1,2);plot(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{2})(locLast1:locLast2));
legend(fields{2},fields{11});title(trials{10})
xlabel('Time (s)');ylabel('Vertical GRF (N)');


% 3. find peak/min, etc of the var so make sure which var you conside for all 4 sets of each trial e.g. S01-04
first = GRFs(ind).data.(fields{2})(locIni1:locIni2);
last = GRFs(ind).data.(fields{2})(locLast1:locLast2);
grfAll =GRFs(ind).data.(fields{2});
% bring any other varaibles here liek SO IK etc.. then do min max peak etc
firstID = IDs(indID).data.knee_angle_l_moment(round(locIni1/20):round(locIni2/20));
lastID = IDs(indID).data.knee_angle_l_moment(round(locLast1/20):round(locLast2/20));
IDhip =  IDs(indID).data.hip_flexion_l_moment;
IDknee =  IDs(indID).data.knee_angle_l_moment;
IDankle =  IDs(indID).data.ankle_angle_l_moment;
IKhip =  IKs(indIK).data.hip_flexion_l;
IKknee =  IKs(indIK).data.knee_angle_l;
IKankle =  IKs(indIK).data.ankle_angle_l;
% fields = fieldnames(SOs(6).data);



% all all muscles etc needed here
[muscleforce] = muslceForces(SOs, string(trials{j}) , 'lquad');
[lquadf(j) ,~] = max(abs(muscleforce));
[muscleforce] = muslceForces(SOs, string(trials{j}) , 'rquad');
[rquadf(j) ,~] = max(abs(muscleforce));
[muscleforce] = muslceForces(SOs, string(trials{j}) , 'lhams');
[lhamsf(j) ,~] = max(abs(muscleforce));
[muscleforce] = muslceForces(SOs, string(trials{j}) , 'rhams');
[rhamsf(j) ,~] = max(abs(muscleforce));

[muscleforce] = muslceForces(SOs, string(trials{j}) , 'lgluts');
[lglutf(j) ,~] = max(abs(muscleforce));
[muscleforce] = muslceForces(SOs, string(trials{j}) , 'rgluts');
[rglutf(j) ,~] = max(abs(muscleforce));

[muscleforce] = muslceForces(SOs, string(trials{j}) , 'lgast');
[lgastf(j) ,~] = max(abs(muscleforce));
[muscleforce] = muslceForces(SOs, string(trials{j}) , 'rgast');
[rgastf(j) ,~] = max(abs(muscleforce));

[muscleforce] = muslceForces(SOs, string(trials{j}) , 'lsole');
[lsolef(j) ,~] = max(abs(muscleforce));
[muscleforce] = muslceForces(SOs, string(trials{j}) , 'rsole');
[rsolef(j) ,~] = max(abs(muscleforce));

[muscleforce] = muslceForces(SOs, string(trials{j}) , 'lTAs');
[lTAsf(j) ,~] = max(abs(muscleforce));
[muscleforce] = muslceForces(SOs, string(trials{j}) , 'rTAs');
[rTAsf(j) ,~] = max(abs(muscleforce));



   % get all gait event vars:
   % lsf = strideFreqL, rsf = strideFreqR, lsl =stepLengthL, rsl:stepLengthR,lsw:stepWidthL, rsw: stepWidthR
   
   %    gaitEvents(indEvent).data.stepWidthR
   [lsf(j)] = mean(gaitEvents(indEvent).data.strideFreqL);%[lsfSD(j)] = std(gaitEvents(indEvent).data.strideFreqL);
   [rsf(j)] = mean(gaitEvents(indEvent).data.strideFreqR);%[rsfSD(j)] = std(gaitEvents(indEvent).data.strideFreqR);
   [lsl(j)] = mean(gaitEvents(indEvent).data.stepLengthL);%[lslSD(j)] = std(gaitEvents(indEvent).data.stepLengthL);
   [rsl(j)] = mean(gaitEvents(indEvent).data.stepLengthR);%[rslSD(j)] = std(gaitEvents(indEvent).data.stepLengthR);
   [lsw(j)] = mean(gaitEvents(indEvent).data.stepWidthL);%[lswSD(j)] = std(gaitEvents(indEvent).data.stepWidthL);
   [rsw(j)] = mean(gaitEvents(indEvent).data.stepWidthR);%[rswSD(j)] = std(gaitEvents(indEvent).data.stepWidthR);

   
   % get all gait Utls vars:
   % bosML = bos_ML, bosAP = bos_AP, mosML =mos_ML, mosAP:mos_AP
%    gaitUtils.bos_ML, gaitUtils.bos_AP, gaitUtils.mos_ML, gaitUtils.mos_AP
   %    gaitEvents(indEvent).data.stepWidthR
   [bosML(j)] = mean(gaitUtils.bos_ML(:,indUtils));%[lsfSD(j)] = std(gaitEvents(indEvent).data.strideFreqL);
   [bosAP(j)] = mean(gaitUtils.bos_AP(:,indUtils));%[rsfSD(j)] = std(gaitEvents(indEvent).data.strideFreqR);
   [mosML(j)] = mean(gaitUtils.mos_ML(:,indUtils));%[lslSD(j)] = std(gaitEvents(indEvent).data.stepLengthL);
   [mosAP(j)] = mean(gaitUtils.mos_AP(:,indUtils));%[rslSD(j)] = std(gaitEvents(indEvent).data.stepLengthR);
   
   
   % get all IK data hip (3dof), knee and ankle
   [hip_flexion_r(j)] = max(abs(IKs(indIK).data.hip_flexion_r));%[lsfSD(j)] = std(gaitEvents(indEvent).data.strideFreqL);
   [hip_adduction_r(j)] = max(abs(IKs(indIK).data.hip_adduction_r));%[lsfSD(j)] = std(gaitEvents(indEvent).data.strideFreqL);
   [hip_rotation_r(j)] =  max(abs(IKs(indIK).data.hip_rotation_r));%[lsfSD(j)] = std(gaitEvents(indEvent).data.strideFreqL);
   [hip_flexion_l(j)] = max(abs(IKs(indIK).data.hip_flexion_l));%[lsfSD(j)] = std(gaitEvents(indEvent).data.strideFreqL);
   [hip_adduction_l(j)] = max(abs(IKs(indIK).data.hip_adduction_l));%[lsfSD(j)] = std(gaitEvents(indEvent).data.strideFreqL);
   [hip_rotation_l(j)] =  max(abs(IKs(indIK).data.hip_rotation_l));%[lsfSD(j)] = std(gaitEvents(indEvent).data.strideFreqL);
   [knee_angle_r(j)] = max(abs(IKs(indIK).data.knee_angle_r));%[lsfSD(j)] = std(gaitEvents(indEvent).data.strideFreqL);
   [knee_angle_l(j)] =  max(abs(IKs(indIK).data.knee_angle_l));
   [ankle_angle_r(j)] = max(abs(IKs(indIK).data.ankle_angle_r));%[lsfSD(j)] = std(gaitEvents(indEvent).data.strideFreqL);
   [ankle_angle_l(j)] =  max(abs(IKs(indIK).data.ankle_angle_l));
   
%    IDhip =  IDs(indID).data.hip_flexion_l_moment;
% IDknee =  IDs(indID).data.knee_angle_l_moment;
% IDankle =  IDs(indID).data.ankle_angle_l_moment;
% IKhip =  IKs(indIK).data.hip_flexion_l;
% IKknee =  IKs(indIK).data.knee_angle_l;
% IKankle =  IKs(indIK).data.ankle_angle_l;
   
   
% Table(1,1) = maxGRF1;Table(1,2) = maxGRF2;

% 4. find mean ± SD or mean(SD) , code this though for a table
% ...

end

% this is the step to be done Tuesda
p=1; % get SD and mean for each grup triasl correct this part
for k=1:4:length(trials)
   lquadMean(p) = mean(lquadf(k:k+3),'omitnan');lquadSD(p) = std(lquadf(k:k+3),'omitnan');
   rquadMean(p) = mean(rquadf(k:k+3),'omitnan');rquadSD(p) = std(rquadf(k:k+3),'omitnan');
   lhamsMean(p) = mean(lhamsf(k:k+3),'omitnan');lhamsSD(p) = std(lhamsf(k:k+3),'omitnan');
   rhamsMean(p) = mean(rhamsf(k:k+3),'omitnan');rhamsSD(p) = std(rhamsf(k:k+3),'omitnan');
   lglutMean(p) = mean(lglutf(k:k+3),'omitnan');llglutSD(p) = std(lglutf(k:k+3),'omitnan');
   rglutMean(p) = mean(rglutf(k:k+3),'omitnan');rlglutSD(p) = std(rglutf(k:k+3),'omitnan');
   
   lgastMean(p) = mean(lgastf(k:k+3),'omitnan');lgastSD(p) = std(lgastf(k:k+3),'omitnan');
   rgastMean(p) = mean(rgastf(k:k+3),'omitnan');rgastSD(p) = std(rgastf(k:k+3),'omitnan');
   lsoleMean(p) = mean(lsolef(k:k+3),'omitnan');lsoleSD(p) = std(lsolef(k:k+3),'omitnan');
   rsoleMean(p) = mean(rsolef(k:k+3),'omitnan');rsoleSD(p) = std(rsolef(k:k+3),'omitnan');
   lTAsMean(p) = mean(lTAsf(k:k+3),'omitnan');lTAsSD(p) = std(lTAsf(k:k+3),'omitnan');
   rTAsMean(p) = mean(rTAsf(k:k+3),'omitnan');rTAsSD(p) = std(rTAsf(k:k+3),'omitnan');
   
     % lsf = strideFreqL, rsf = strideFreqR, lsl =stepLengthL, rsl:stepLengthR,lsw:stepWidthL, rsw: stepWidthR
   
   %    gaitEvents(indEvent).data.stepWidthR
   [lsfMean(p)] = mean(lsf(k:k+3),'omitnan');[lsfSD(p)] = std(lsf(k:k+3),'omitnan');
   [rsfMean(p)] = mean(rsf(k:k+3),'omitnan');[rsfSD(p)] = std(lsf(k:k+3),'omitnan');
   [lslMean(p)] = mean(lsl(k:k+3),'omitnan');[lslSD(p)] = std(lsf(k:k+3),'omitnan');
   [rslMean(p)] = mean(rsl(k:k+3),'omitnan');[rslSD(p)] = std(lsf(k:k+3),'omitnan');
   [lswMean(p)] = mean(lsw(k:k+3),'omitnan');[lswSD(p)] = std(lsf(k:k+3),'omitnan');
   [rswMean(p)] = mean(rsw(k:k+3),'omitnan');[rswSD(p)] = std(lsf(k:k+3),'omitnan');
     
   
   % add the MoS, Com, bos, etc here:
%    bosML = bos_ML, bosAP = bos_AP, mosML =mos_ML, mosAP:mos_AP
   [bosMLMean(p)] = mean(bosML(k:k+3),'omitnan');[bosMLSD(p)] = std(bosML(k:k+3),'omitnan');
   [bosAPMean(p)] = mean(bosAP(k:k+3),'omitnan');[bosAPSD(p)] = std(bosAP(k:k+3),'omitnan');
   [mosMLMean(p)] = mean(mosML(k:k+3),'omitnan');[mosMLSD(p)] = std(mosML(k:k+3),'omitnan');
   [mosAPMean(p)] = mean(mosAP(k:k+3),'omitnan');[mosAPSD(p)] = std(mosAP(k:k+3),'omitnan');
      
   
   % IK outcomes
   
   [hip_flexion_rMean(p)] = mean(hip_flexion_r(k:k+3),'omitnan');[hip_flexion_rSD(p)] = std(hip_flexion_r(k:k+3),'omitnan');
   [hip_flexion_lMean(p)] = mean(hip_flexion_l(k:k+3),'omitnan');[hip_flexion_lSD(p)] = std(hip_flexion_l(k:k+3),'omitnan');
   [hip_adduction_rMean(p)] = mean(hip_adduction_r(k:k+3),'omitnan');[hip_adduction_rSD(p)] = std(hip_adduction_r(k:k+3),'omitnan');
   [hip_adduction_lMean(p)] = mean(hip_adduction_l(k:k+3),'omitnan');[hip_adduction_lSD(p)] = std(hip_adduction_l(k:k+3),'omitnan');
   [hip_rotation_rMean(p)] = mean(hip_rotation_r(k:k+3),'omitnan');[hip_rotation_rSD(p)] = std(hip_rotation_r(k:k+3),'omitnan');
   [hip_rotation_lMean(p)] = mean(hip_rotation_l(k:k+3),'omitnan');[hip_rotation_lSD(p)] = std(hip_rotation_l(k:k+3),'omitnan');
   [knee_angle_rMean(p)] = mean(knee_angle_r(k:k+3),'omitnan');[knee_angle_rSD(p)] = std(knee_angle_r(k:k+3),'omitnan');
   [knee_angle_lMean(p)] = mean(knee_angle_l(k:k+3),'omitnan');[knee_angle_lSD(p)] = std(knee_angle_l(k:k+3),'omitnan');
   [ankle_angle_rMean(p)] = mean(ankle_angle_r(k:k+3),'omitnan');[ankle_angle_rSD(p)] = std(ankle_angle_r(k:k+3),'omitnan');
   [ankle_angle_lMean(p)] = mean(ankle_angle_l(k:k+3),'omitnan');[ankle_angle_lSD(p)] = std(ankle_angle_l(k:k+3),'omitnan');
   
   p=p+1;
end

% 5. Put it in a table format... create table 
Variables = {'Left Quad','Right Quad','Left Hams','Right Hams',...
    'Left Glut','Right Glut','Left Gast','Right Gast','Left Sol','Right Sol','Left TA','Right TA',...
    'Left Step Frequency','Right Step Frequency','Left Step Length','Right Step Length',...
    'Left Step Width','Right Step Width','MoS AP','MoS ML', 'BoS AP','BoS ML'...
    ,'Right Hip Flexion','Left Hip Flexion','Right Hip Adduction','Left Hip Adduction',...
    'Right Hip Rotation','Left Hip Rotation','Right Knee Flexion','Left Knee Flexion', ...
    'Right Ankle Flexion','Left Ankle Flexion'};
%     Peak vGRF';'Hip Joint Moment';'Knee Joint Moment';'Ankle Joint Moment';...
%     'Hamstrings';'Quadriceps'};
% ;'MoS A/P';'MoS M/L'; 'Step Length';'Step Frequency';...
%    'Cost of Dual Tasking';'Cost';'C2ost';'IDall';'SDIDall'}; % rows
ColName =  {'S';'W1';'W1+T';'W2';'W2+T';'W3';'W3+T'};
S = [lquadMean(1);lquadSD(1);lquadMean(1);lquadSD(1);lquadMean(1);lquadSD(1);0;0;0;0;0;0];
W1 =S ;
W1T =S;
W2 = S;
W2T = S;
W3 = S;
W3T = S;
% T = table(S,W1,W1T,W2,W2T,W3,W3T,'RowNames',Variables);
% uitable('Data',T{:,:},'ColumnName',ColName,...
%      'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

% 6. save the table and save it as figure for the final publication
% chec this below:
% https://www.mathworks.com/matlabcentral/answers/254690-how-can-i-display-a-matlab-table-in-a-figure
%  mean/sd in a table, see thisL
% https://www.mathworks.com/matlabcentral/answers/52731-display-mean-std-in-a-table-or-matrix

% this must be the last step when we have all the mean and SDs of all
% varaibes
% Create sample data and row and column headers.

hh= {lquadMean; rquadMean;lhamsMean;rhamsMean;lglutMean;rglutMean;lgastMean;...
    rgastMean;lsoleMean;rsoleMean;lTAsMean;rTAsMean;lsfMean;rsfMean;lslMean;rslMean...
    ;lswMean;rswMean;mosAPMean;mosMLMean;bosAPMean;bosMLMean;hip_flexion_rMean...
    ;hip_flexion_lMean;hip_adduction_rMean;hip_adduction_lMean;hip_rotation_rMean;...
    hip_rotation_lMean;knee_angle_rMean;knee_angle_lMean;ankle_angle_rMean;ankle_angle_lMean};


pp={lquadSD;rquadSD;lhamsSD;rhamsSD;llglutSD;rlglutSD;lgastSD;rgastSD;lsoleSD...
    ;rsoleSD;lTAsSD;rTAsSD;lsfSD;rsfSD;lslSD;rslSD;lswSD;rswSD;mosAPSD;mosMLSD;bosAPSD;bosMLSD;...
    hip_flexion_rSD;hip_flexion_lSD;hip_adduction_rSD;hip_adduction_lSD;hip_rotation_rSD;...
    hip_rotation_lSD;knee_angle_rSD;knee_angle_lSD;ankle_angle_rSD;ankle_angle_rSD};
 

figure
columnHeaders = {'n', 'Data Set #1', 'Data Set #2'};
columnHeaders = ColName;
for t=1:length(hh)
for n=1:length(trials)/4
%   rowHeaders{n} = sprintf('Row #%d', n);
%   tableData{n,1} = n;
  tableData{t,n} = sprintf('%.2f %s %.2f', hh{t}(n), 177, pp{t}(n));%10*rand(1,1);
%   tableData{3,n} = sprintf('%.2f %s %.2f', hh{t}(n), 177, pp{t}(n));
end
end
% Create the table and display it.
hTable = uitable();
% Apply the row and column headers.
rowHeaders = Variables;
set(hTable, 'RowName', rowHeaders);
set(hTable, 'ColumnName', columnHeaders);
% Display the table of values.
set(hTable, 'data', tableData);
% Size the table.
set(hTable, 'units', 'normalized');
set(hTable, 'Position', [.1 .1 .8 .8]);
set(hTable, 'Position', [.1 .1 .85 .85]);
set(hTable, 'ColumnWidth', {40, 120, 180});
% set(hTable, 'ColumnWidth', {80,70,70,70});
w = 90;
hTable.ColumnWidth{1} = w;
hTable.ColumnWidth{2} = w;
hTable.ColumnWidth{3} = w;
hTable.ColumnWidth{4} = w;
hTable.ColumnWidth{5} = w;
hTable.ColumnWidth{6} = w;
hTable.ColumnWidth{7} = w;
set(gcf,'name','Image Analysis Demo','numbertitle','off')



%% cost of dual tasking + perforance 4 graphs
% (1,4): including hit performance + cost of it
% overal performance + cost of it
source =  'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\Publications\Biomecognition_PaperMethod\Results';
load(fullfile(source,'SH.mat')) % import SH all data for performance etc
%
% lets load the SH data from experiemnt and trigger etc and then plot the performance 
% load('C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\Publications\Biomecognition_PaperMethod\Results\SH')
clear a b
subplot(2,2,1)
k=1;
% for i = 1:length(SH.data.TrialNo)
for i = [2, 3, 4, 7, 8, 11, 12] % for the paper just dual task:[2 3 4 7 8 11 12]
if ~isnan(SH.data.Performance(i))
    a(k)= SH.data.Performance(i);
    b(k)= SH.data.TrialType(i);
k=k+1;
end
end
c=categorical(string(b));
bar(c,a)
xlabel('Trials');ylabel('Overall Performance (%)')
ylim([0 100])

% cost calucaltion with respect to T only
subplot(2,2,3)
cost1 = (a(1)-a)/a(1)*100;
c=categorical(string(b));
bar(c,cost1)
xlabel('Trials');ylabel('Cost of Overall Performance (%)')
ylim([-20 100])

clear a b
subplot(2,2,2)
k=1;
% for i = 1:length(SH.data.TrialNo)
for i = [2, 3, 4, 7, 8, 11, 12] % for the paper just dual task:[2 3 4 7 8 11 12]
if ~isnan(SH.data.HitPerformance(i))
    a(k)= SH.data.HitPerformance(i);
    b(k)= SH.data.TrialType(i);
k=k+1;
end
end
c=categorical(string(b));
bar(c,a)
xlabel('Trials');ylabel('Hit Performance (%)')
ylim([0 100])

% cost calucaltion with respect to T only
subplot(2,2,4)
cost1 = (a(1)-a)/a(1)*100;
c=categorical(string(b));
bar(c,cost1)
xlabel('Trials');ylabel('Cost of Hit Performance (%)')
ylim([-20 100])


%% GRF all 7s and averaged for 4 intervals of 7s.
% I don't think this make sence we need to have FS to FS and first and last would be perfect so ignore the idea of 7s

% find initial, end pertub from dflow, then find two gait cycles: (HS-TO ... Preturb ...HS-TO)
clear 
close all
clc

vars = {'vGRF'}; % then add any variable from the table

% 0. Path and data location etc and whic trials

addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';
load(fullfile(source, 'DFLOW\dflow.mat'))
load(fullfile(source, 'trcResults\gaitEvents.mat'))
trials = {'S+T01','S+T02','S+T03','S+T04', ...
'W101','W102','W103','W104',...
    'W1+T01','W1+T02','W1+T03','W1+T04',...
    'W201','W202','W203','W204',...
    'W2+T01','W2+T02','W2+T03','W2+T04',...
    'W301','W302','W303','W304',...
    'W3+T01','W3+T02','W3+T03','W3+T04'};

% load all the data GRF, muscles, mos, etc
load(fullfile(source,'grfResults\GRFs.mat'))
% load(fullfile(source,'SOResults\SOs.mat'))
% load(fullfile(source,'IDResults\IDs.mat'))
% load(fullfile(source,'IKResults\IKs.mat'))

k=1;
for j=1:length(trials)
% 1. load vars{1} e.g. vGRF load GRFs

load(fullfile(source,'grfResults\GRFs.mat'))

fields = fieldnames(GRFs(6).data);
for i = 1:length(GRFs)
    trials2{i}= GRFs(i).name;
%     try 
%     trialsID{i}= extractBefore(IDs(i).name,'_');
%     trialsSO{i}=  extractBefore(SOs(i).name,'_');
%     trialsIK{i}= extractBefore(IKs(i).name,'_');
%     catch 
%         fprintf('Inconsistent data in iteration %s, skipped.\n', i);
%     end
end

[ind(k), missTrials]=trialIndexFinder(string(trials{j}),trials2);
k=k+1;
% try 
% [indSO, missTrialsSO]=trialIndexFinder(string(trials{j}),trialsSO);
% [indID, missTrialsID]=trialIndexFinder(string(trials{j}),trialsID);
% [indIK, missTrialsIK]=trialIndexFinder(string(trials{j}),trialsIK);
% catch 
%     print('chcek later')
% end
% ini = dflow.Perturb(28).iniPerturb;

end


% plot the average of 4 sets of 7 second in every 2 min.
p=1; % get SD and mean for each grup triasl
for k=1:4:length(trials)
figure
       a1 = (GRFs(ind(k)).data.(fields{2}));
   a2 = (GRFs(ind(k+1)).data.(fields{2}));
   a3 = (GRFs(ind(k+2)).data.(fields{2}));
   a4 = (GRFs(ind(k+3)).data.(fields{2}));
%    
hold on;plot(a1);title(GRFs(ind(k)).name)
plot(a2);title(GRFs(ind(k+1)).name)
plot(a3);title(GRFs(ind(k+2)).name)
plot(a4);title(GRFs(ind(k+3)).name)
%    s = max([size(a1) size(a2) size(a3) size(a4)]);
%    a1 = ZeroTo100(a1,[],s);a2 = ZeroTo100(a2,[],s);a3 = ZeroTo100(a3,[],s);a4 = ZeroTo100(a4,[],s);
%   meanGrf(:,p) = mean([a1,a2,a3,a4],2);
%    p=p+1;
end

h = figure;
x = cumsum(randn(10,95),2);
y = cumsum(randn(10,95),2)*2;
z = cumsum(rand(10,95),2)/2;
stdshade(x,0.5,'b');
hold(h.Children,'on');
stdshade(y,0.5,'r');
stdshade(z,0.5,'g');
axis square; xlim([1 95]);

subplot(2,1,1);plot(GRFs(ind).data.(fields{end}),GRFs(ind).data.(fields{2}));
subplot(2,1,2);plot(GRFs(ind).data.(fields{end}),GRFs(ind).data.(fields{2}));


legend(fields{2},fields{11});title(trials{10})

xlabel('Time (s)');ylabel('Vertical GRF (N)')


%% •GRF
% an example of GRF with TO for left side 
figure
import org.opensim.modeling.*
cd(fullfile(source,'grfResults'))
timeSeriesTable = STOFileAdapter.read(trialgrf); 
% STOFileAdapter.write(timeSeriesTable,'test.mot');


Forces = osimTableToStruct(timeSeriesTable);


fields = fieldnames(Forces);
% check if there is nan in the mot files
% fnames = findnanfields(Forces);


% figure
% for i = 1:length(fields)-1
% hold on
% 
%     plot(Forces.(fields{i}))
% end
% figure
 hold on;plot(Forces.(fields{end}),Forces.(fields{2}));plot(Forces.(fields{end}),Forces.(fields{11}))
 legend(fields{2},fields{11});title(trialgrf)
 xlabel('Time (s)');ylabel('Vertical GRF (N)')

%% •Joint Kinematics and kinetics
% IK reslts plot
figure
import org.opensim.modeling.*
cd(fullfile(source,'IKResults'))
IKTable = STOFileAdapter.read(trialIK);
IKStruct = osimTableToStruct(IKTable); 


figure
fields = fieldnames(IKStruct);
for i = 1:length(fields)-1
hold on
    plot(IKStruct.(fields{i}))
end

figure
hold on; plot(IKStruct.(fields{10}));plot(IKStruct.(fields{17}))
legend(fields{10},fields{17})
figure
hold on; plot(IKStruct.(fields{7}));plot(IKStruct.(fields{14}))
legend(fields{7},fields{14})

%% moments
%% ID reslts plot % GRF dont seem to be wokring check this well
figure
import org.opensim.modeling.*
cd(fullfile(source,'IDResults'))
IDTable = STOFileAdapter.read(trialID);
IDStruct = osimTableToStruct(IDTable);
fields = fieldnames(IDStruct);


for i = 1:length(fields)-1
hold on
    plot(IDStruct.(fields{i}))
end
figure
hold on; plot(-IDStruct.knee_angle_r_moment);plot(-IDStruct.knee_angle_l_moment)
legend('Right Knee Moment','Left Knee Moment')


%% •Muscle contribution

% SO reslts plot
figure
import org.opensim.modeling.*
cd(fullfile(source,'SOResults'))
SOTable = STOFileAdapter.read(trialSO);
SOStruct = osimTableToStruct(SOTable);
fields = fieldnames(SOStruct);

for i = 1:length(fields)-1
hold on
    plot(SOStruct.(fields{i}))
end
figure
hold on; plot(SOStruct.(fields{28}));plot(SOStruct.(fields{35}))
legend(fields{28},fields{35})

% [8-11] is semimem_r	semiten_r	bifemlh_r	bifemsh_r or right
% Hamstring and left: [51-54]: Right Quad: [29:32], Left Quad: [72:75], 

HamR = 8:11; HamL = 51:54;
QuadR = 29:32; QuadL = 72:75;
RHam = 0; LHam=0; RQuad = 0; LQuad=0;
for i = 1:length(HamR)
Buff = SOStruct.(fields{HamR(i)});
BuffL = SOStruct.(fields{HamL(i)});
RHam = Buff + RHam;
LHam = BuffL + LHam;

BuffQ = SOStruct.(fields{QuadR(i)}); %Right Quad force
BuffLQ = SOStruct.(fields{QuadL(i)}); %Left Quad force
RQuad = BuffQ + RQuad;
LQuad = BuffLQ + LQuad;

end

figure()
hold on; plot(RHam);plot(LHam);
legend('Right Ham Force','Left Ham Force')
title(trialSO)
figure()
hold on; plot(RQuad);plot(LQuad);
legend('Right Quad Force','Left Quad Force')

title(trialSO)

% import org.opensim.modeling.*
% analyzeTool = AnalyzeTool('SO_HBM_Setup.xml');



%% GRF pre post- perturb gait cycle:

% find initial, end pertub from dflow, then find two gait cycles: (HS-TO ... Preturb ...HS-TO)

addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';
load(fullfile(source, 'DFLOW\dflow.mat'))
load(fullfile(source, 'trcResults\gaitEvents.mat'))


figure
import org.opensim.modeling.*
load(fullfile(source,'grfResults\GRFs.mat'))

fields = fieldnames(GRFs(6).data);
for i = 1:length(GRFs)
    trials2{i}= GRFs(i).name;
end

[ind, missTrials]=trialIndexFinder(trialgrf,trials2);
% ini = dflow.Perturb(28).iniPerturb;
% last = dflow.Perturb(28).endPertub;
%lets find first and last cycles
ini1 = gaitEvents(ind).data.HSLocL(1)
ini2 = gaitEvents(ind).data.HSLocL(2)
last1 = gaitEvents(ind).data.HSLocL(end-1)
last2 = gaitEvents(ind).data.HSLocL(end)
[~,locIni1]=min(abs(ini1-GRFs(ind).data.(fields{end})));
[~,locIni2]=min(abs(ini2-GRFs(ind).data.(fields{end})));
[~,locLast1]=min(abs(last1-GRFs(ind).data.(fields{end})));
[~,locLast2]=min(abs(last2-GRFs(ind).data.(fields{end})));
subplot(2,1,1);plot(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{2})(locIni1:locIni2));
subplot(2,1,2);plot(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{2})(locLast1:locLast2));
 legend(fields{2},fields{11});title(trialgrf)
 xlabel('Time (s)');ylabel('Vertical GRF (N)')


 %% GRF first FS to second FS in 7s for all three speeds+T
 
 % let;s plot W101, W201,W301 vs.W1+T01, W2+T01,W3+T01 Left side in Vertcal


% find initial, end pertub from dflow, then find two gait cycles: (HS-TO ... Preturb ...HS-TO)

addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';
load(fullfile(source, 'DFLOW\dflow.mat'))
load(fullfile(source, 'trcResults\gaitEvents.mat'))

figure
import org.opensim.modeling.*
load(fullfile(source,'grfResults\GRFs.mat'))
% grfs = string({'W101','W202','W301.mot','W1+T01.mot', 'W2+T01.mot','W3+T01.mot'});
grfs = string({'W101','W202','W301','W1+T01', 'W2+T01','W3+T01'});

fields = fieldnames(GRFs(6).data); % just to get the name of feilds so no change
% get all the names wihtin GRFs all trials in fact
for i = 1:length(GRFs)
    trials2{i}= GRFs(i).name;
end

%find the index of the ones (grfs) we aim to plot within trails2
[ind1, missTrials]=trialIndexFinder(grfs,trials2);
% ini = dflow.Perturb(28).iniPerturb;
% last = dflow.Perturb(28).endPertub;
%lets find first and last cycles
for ind=ind1%[ind1(3) ind1(end)]%

    ini1 = gaitEvents(ind).data.HSLocL(1)
ini2 = gaitEvents(ind).data.HSLocL(2)
last1 = gaitEvents(ind).data.HSLocL(end-1)
last2 = gaitEvents(ind).data.HSLocL(end)
[~,locIni1]=min(abs(ini1-GRFs(ind).data.(fields{end})));
[~,locIni2]=min(abs(ini2-GRFs(ind).data.(fields{end})));
[~,locLast1]=min(abs(last1-GRFs(ind).data.(fields{end})));
[~,locLast2]=min(abs(last2-GRFs(ind).data.(fields{end})));
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{2})(locIni1:locIni2));
t=[1:1:101];subplot(1,2,1);hold on;plot(t,R);
%  legend(fields{2},fields{11});title(trialgrf)
legend(grfs)
xlabel('Time (s)');ylabel('Vertical GRF (N)')
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{2})(locLast1:locLast2))
xlim([0 100]);title('GRF comparison left side - first gait cycle');
subplot(1,2,2);hold on;plot(t,R);
%  legend(fields{2},fields{11});title(trialgrf)
 xlabel('Time (s)');ylabel('Vertical GRF (N)')
legend(grfs);title('GRF comparison left side - Last gait cycle');
xlim([0 100])
    
end

%% GRF: Cost of dual task = (ST_GRFmax - DT_GRFmax) / ST_GRFmax

% lets just find W101 W1+T01 S+T01 compare them cnside the above and below eq.
% Cost of dual task = (ST_GRFmax - DT_GRFmax) / ST_GRFmax

addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';
load(fullfile(source, 'DFLOW\dflow.mat'))
load(fullfile(source, 'trcResults\gaitEvents.mat'))

% figure
import org.opensim.modeling.*
load(fullfile(source,'grfResults\GRFs.mat'))
grfs = string({'W101.mot','W202.mot','W301.mot','W1+T01.mot', 'W2+T01.mot','W3+T01.mot'});

fields = fieldnames(GRFs(6).data); % just to get the name of feilds so no change
% get all the names wihtin GRFs all trials in fact
for i = 1:length(GRFs)
    trials2{i}= GRFs(i).name;
end

%find the index of the ones (grfs) we aim to plot within trails2
[ind1, missTrials]=trialIndexFinder(grfs,trials2);

k=1;
for ind=ind1%[ind1(3) ind1(end)]%

  peaksGRF(k) = max(GRFs(ind).data.(fields{2}))
  k=k+1;
end
k=1;
for i=1:length(grfs)/2
dtCost(k) = (peaksGRF(i)-peaksGRF(i+3))/peaksGRF(i)*100
k=k+1;
end
bar(c,dtCost)
c=categorical(grfs(1:3))
ylabel('Cost of Dual Task in peak GRF')
xlabel('Different Gait Speeds ')

% bar(peaksGRF)

%% joint angles hip-knee-ankle 0:100% 
% 
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
trials = string({'W101_ik','W202_ik','W301_ik','W1+T01_ik', 'W2+T01_ik','W3+T01_ik'});
trials = string({'W101_ik','W1+T01_ik','W1+T+L+DEC_ik','W1+T+L+ACC_ik','W1+T+R+DEC_ik','W1+T+R+ACC_ik'});
trials = string({'W201_ik','W2+T01_ik','W2+T+L+DEC_ik','W2+T+L+ACC_ik','W2+T+R+DEC_ik','W2+T+R+ACC_ik'});
trials = string({'W301_ik','W3+T01_ik','W3+T+L+DEC_ik','W3+T+L+ACC_ik','W3+T+R+DEC_ik','W3+T+R+ACC_ik'});
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';
load(fullfile(source, 'DFLOW\dflow.mat'))
load(fullfile(source, 'trcResults\gaitEvents.mat'))

load(fullfile(source, 'IKResults\IKs.mat'))
xx = [0 100];
yy = [-100 50];

fields = fieldnames(IKs(6).data);
% for i = 1:length(fields)-1
% hold on
%     plot(IKs.data.(fields{i}))
% end

for i = 1:length(IKs)
    trials2{i}= IKs(i).name;
end

%find the index of the ones (grfs) we aim to plot within trails2
[ind1, missTrials]=trialIndexFinder(trials,trials2);

clear GRFs
GRFs = IKs;

for ind=ind1%[ind1(3) ind1(end)]%

    ini1 = gaitEvents(ind).data.HSLocL(1)
ini2 = gaitEvents(ind).data.HSLocL(2)
last1 = gaitEvents(ind).data.HSLocL(end-1)
last2 = gaitEvents(ind).data.HSLocL(end)
[~,locIni1]=min(abs(ini1-GRFs(ind).data.(fields{end})));
[~,locIni2]=min(abs(ini2-GRFs(ind).data.(fields{end})));
[~,locLast1]=min(abs(last1-GRFs(ind).data.(fields{end})));
[~,locLast2]=min(abs(last2-GRFs(ind).data.(fields{end})));

%get first stance phase
%hip right
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{14})(locIni1:locIni2),100);
% all 7second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{2}));
t=[1:1:101];subplot(3,2,1);hold on;plot(t,R);title('First Stance Phase (%)');
xlim(xx)
ylim(yy)
%  legend(fields{2},fields{11});title(trialgrf)
% legend(trials)
xlabel('Time (s)');ylabel(' Left Hip Flexion (deg)')

%knee right
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{17})(locIni1:locIni2),100);
% all 7second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{2}));
t=[1:1:101];subplot(3,2,3);hold on;plot(t,R);
xlim(xx)
ylim(yy)
%  legend(fields{2},fields{11});title(trialgrf)
 legend(trials)
xlabel('Time (s)');ylabel(' Left Knee Flexion (deg)')

%ankle right
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{18})(locIni1:locIni2),100);
% all 7second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{2}));
t=[1:1:101];subplot(3,2,5);hold on;plot(t,R);
xlim(xx)
ylim(yy)
%  legend(fields{2},fields{11});title(trialgrf)
% legend(trials)
xlabel('Time (s)');ylabel(' Left Ankle Flexion (deg)')





%get last stance phase
%hip left
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{14})(locLast1:locLast2),100)
%all 7 second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{7}))

subplot(3,2,2);hold on;plot(t,R);title('Last Stance Phase (%)');
%  legend(fields{2},fields{11});title(trialgrf)
 xlabel('Time (s)');ylabel(' Left Hip Flexion (deg)')
% legend(trials)
xlim(xx)
ylim(yy)
    

%Knee left
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{17})(locLast1:locLast2),100)
%all 7 second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{7}))
xlim(xx)
ylim(yy)
subplot(3,2,4);hold on;plot(t,R);
%  legend(fields{2},fields{11});title(trialgrf)
 xlabel('Time (s)');ylabel(' Left Knee Flexion (deg)')
% legend(trials)
xlim(xx)
ylim(yy)


%hip left
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{18})(locLast1:locLast2),100)
%all 7 second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{7}))
xlim(xx)
ylim(yy)
subplot(3,2,6);hold on;plot(t,R);
%  legend(fields{2},fields{11});title(trialgrf)
 xlabel('Time (s)');ylabel(' Left Ankle Flexion (deg)')
% legend(trials)
xlim(xx)
ylim(yy)


end





% 
% figure
% hold on; plot(IKs.(fields{10}));plot(IKs.(fields{17}))
% legend(fields{10},fields{17})
% figure
% hold on; plot(IKs.(fields{7}));plot(IKs.(fields{14}))
% legend(fields{7},fields{14})


%% joint moments angles hip-knee-ankle 0:100% 
% 
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
trials = string({'W101_ID','W202_ID','W302_ID','W1+T01_ID', 'W2+T01_ID','W3+T01_ID'});
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';
load(fullfile(source, 'DFLOW\dflow.mat'))
load(fullfile(source, 'trcResults\gaitEvents.mat'))

load(fullfile(source, 'IDResults\IDs.mat'))
xx = [0 100];
yy = [-140 120];

fields = fieldnames(IDs(6).data);
% for i = 1:length(fields)-1
% hold on
%     plot(IKs.data.(fields{i}))
% end

for i = 1:length(IDs)
    trials2{i}= IDs(i).name;
end

%find the index of the ones (grfs) we aim to plot within trails2
[ind1, missTrials]=trialIndexFinder(trials,trials2);
trials = string({'W101\_ID','W202\_ID','W302\_ID','W1+T01\_ID', 'W2+T01\_ID','W3+T01\_ID'});

clear GRFs
GRFs = IDs;

for ind=ind1%[ind1(3) ind1(end)]%

ini1 = gaitEvents(ind).data.HSLocL(1)
ini2 = gaitEvents(ind).data.HSLocL(2)
last1 = gaitEvents(ind).data.HSLocL(end-1)
last2 = gaitEvents(ind).data.HSLocL(end)
[~,locIni1]=min(abs(ini1-GRFs(ind).data.(fields{end})));
[~,locIni2]=min(abs(ini2-GRFs(ind).data.(fields{end})));
[~,locLast1]=min(abs(last1-GRFs(ind).data.(fields{end})));
[~,locLast2]=min(abs(last2-GRFs(ind).data.(fields{end})));

%get first stance phase
%hip right
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{10})(locIni1:locIni2),100);
% all 7second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{2}));
t=[1:1:101];subplot(3,2,1);hold on;plot(t,R);title('First Stance Phase (%)');
xlim(xx)
ylim(yy)
%  legend(fields{2},fields{11});title(trialgrf)
% legend(trials)
xlabel('Time (s)');ylabel(' Left Hip Flexion (N.m)')

%knee right
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{17})(locIni1:locIni2),100);
% all 7second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{2}));
t=[1:1:101];subplot(3,2,3);hold on;plot(t,R);
xlim(xx)
ylim(yy)
%  legend(fields{2},fields{11});title(trialgrf)
 legend(trials)
xlabel('Time (s)');ylabel(' Left Knee Flexion (N.m)')

%ankle right
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{19})(locIni1:locIni2),100);
% all 7second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{2}));
t=[1:1:101];subplot(3,2,5);hold on;plot(t,R);
xlim(xx)
ylim(yy)
%  legend(fields{2},fields{11});title(trialgrf)
% legend(trials)
xlabel('Time (s)');ylabel(' Left Ankle Flexion (N.m)')



%get last stance phase
%hip left
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{10})(locLast1:locLast2),100)
%all 7 second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{7}))

subplot(3,2,2);hold on;plot(t,R);title('Last Stance Phase (%)');
%  legend(fields{2},fields{11});title(trialgrf)
 xlabel('Time (s)');ylabel(' Left Hip Flexion (N.m)')
% legend(trials)
xlim(xx)
ylim(yy)
    

%Knee left
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{17})(locLast1:locLast2),100)
%all 7 second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{7}))
xlim(xx)
ylim(yy)
subplot(3,2,4);hold on;plot(t,R);
%  legend(fields{2},fields{11});title(trialgrf)
 xlabel('Time (s)');ylabel(' Left Knee Flexion (N.m)')
% legend(trials)
xlim(xx)
ylim(yy)


%hip left
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{19})(locLast1:locLast2),100)
%all 7 second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{7}))
xlim(xx)
ylim(yy)
subplot(3,2,6);hold on;plot(t,R);
%  legend(fields{2},fields{11});title(trialgrf)
 xlabel('Time (s)');ylabel(' Left Ankle Flexion (N.m)')
% legend(trials)
xlim(xx)
ylim(yy)


end



%% Muscle quad and ham0:100% first and last cycle
% Lets do these ones: Quad, Ham, Sol, Gas, Glut, TA 
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
trials = string({'W101_StaticOptimization_force','W201_StaticOptimization_force','W302_StaticOptimization_force',...
    'W1+T01_StaticOptimization_force', 'W2+T01_StaticOptimization_force','W3+T01_StaticOptimization_force'});
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\CB1\';
load(fullfile(source, 'DFLOW\dflow.mat'))
load(fullfile(source, 'trcResults\gaitEvents.mat'))

load(fullfile(source, 'SOResults\SOs.mat'))
xx = [0 100];
yy = [-0 2000];

fields = fieldnames(SOs(6).data);
% for i = 1:length(fields)-1
% hold on
%     plot(IKs.data.(fields{i}))
% end

for i = 1:length(SOs)
    trials2{i}= SOs(i).name;
end

%find the index of the ones (grfs) we aim to plot within trails2
[ind1, missTrials]=trialIndexFinder(trials,trials2);
trials1 = string({'W101\_SO','W201\_SO','W302\_SO','W1+T01\_SO', 'W2+T01\_SO','W3+T01\_SO'});

clear GRFs
GRFs = SOs;


HamR = 7:10; HamL = 51:54;
QuadR = 29:32; QuadL = 72:75;
SolR = 35; SolL = 78;
GasR = 33:34; GasL = 76:77;
GlutR = [2:7 21:23]; GlutL = [45:50 64:66];
TAR = 39; TAL = 82;


RHam = 0; LHam=0; RQuad = 0; LQuad=0;
RSol = 0; LSol=0; RGas = 0; LGas=0;
RGlut = 0; LGlut=0; RTA = 0; LTA=0;



for ind=ind1%[ind1(3) ind1(end)]%

ini1 = gaitEvents(ind).data.HSLocL(1)
ini2 = gaitEvents(ind).data.HSLocL(2)
last1 = gaitEvents(ind).data.HSLocL(end-1)
last2 = gaitEvents(ind).data.HSLocL(end)
[~,locIni1]=min(abs(ini1-GRFs(ind).data.(fields{end})));
[~,locIni2]=min(abs(ini2-GRFs(ind).data.(fields{end})));
[~,locLast1]=min(abs(last1-GRFs(ind).data.(fields{end})));
[~,locLast2]=min(abs(last2-GRFs(ind).data.(fields{end})));

% this was to get Quad and Ham in R/L for those we want to plot
for i = 1:length(HamR)
Buff = SOs(ind).data.(fields{HamR(i)});
BuffL = SOs(ind).data.(fields{HamL(i)});
RHam = Buff + RHam;
LHam = BuffL + LHam;

BuffQ = SOs(ind).data.(fields{QuadR(i)}); %Right Quad force
BuffLQ = SOs(ind).data.(fields{QuadL(i)}); %Left Quad force
RQuad = BuffQ + RQuad;
LQuad = BuffLQ + LQuad;
end
% Sol, 
BuffS = SOs(ind).data.(fields{SolR(1)}); %Right Quad force
BuffLS = SOs(ind).data.(fields{SolL(1)}); %Left Quad force
RSol = BuffS + RSol;
LSol = BuffLS + LSol;

%  Gas, 
for i = 1:length(GasR)
BuffGa = SOs(ind).data.(fields{GasR(i)}); %Right Quad force
BuffLGa = SOs(ind).data.(fields{GasL(i)}); %Left Quad force
RGas = BuffGa + RGas;
LGas = BuffLGa + LGas;
end

%  Glut,  
for i = 1:length(GlutR)
BuffG = SOs(ind).data.(fields{GlutR(i)}); %Right Quad force
BuffLG = SOs(ind).data.(fields{GlutL(i)}); %Left Quad force
RGlut = BuffG + RGlut;
LGlut = BuffLG + LGlut;
end
%  TA 
BuffA = SOs(ind).data.(fields{TAR(1)}); %Right Quad force
BuffLA = SOs(ind).data.(fields{TAL(1)}); %Left Quad force
RTA = BuffA + RTA;
LTA = BuffLA + LTA;




%get first stance phase
%hip right
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),LHam(locIni1:locIni2),100);
% all 7second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{2}));
t=[1:1:101];subplot(3,2,1);hold on;plot(t,R);title('Ham muslce Left, First Gait Cycle (%)')
xlim(xx)
ylim(yy)
%  legend(fields{2},fields{11});title(trialgrf)
% legend(trials)
xlabel('Time (s)');ylabel(' Muscle force (N)')

%knee right
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),LQuad(locIni1:locIni2),100);
% all 7second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{2}));
t=[1:1:101];subplot(3,2,3);hold on;plot(t,R);title('Quad muslce Left');
xlim(xx)
ylim(yy)
%  legend(fields{2},fields{11});title(trialgrf)
 
xlabel('Time (s)');ylabel(' Muscle force (N)')

%ankle right
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),LSol(locIni1:locIni2),100);
% all 7second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locIni1:locIni2),GRFs(ind).data.(fields{2}));
t=[1:1:101];subplot(3,2,5);hold on;plot(t,R);title('Sol muslce Left');
xlim(xx)
ylim(yy)
%  legend(fields{2},fields{11});title(trialgrf)
% legend(trials)
xlabel('Time (s)');ylabel(' Muscle force (N)')
legend(trials1)


%get last stance phase
%hip left
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),LHam(locLast1:locLast2),100)
%all 7 second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{7}))

subplot(3,2,2);hold on;plot(t,R);title('Ham muslce Left, Last Gait Cycle (%)');
%  legend(fields{2},fields{11});title(trialgrf)
 xlabel('Time (s)');ylabel(' Muscle force (N)')
% legend(trials)
xlim(xx)
ylim(yy)
    

%Knee left
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),LQuad(locLast1:locLast2),100)
%all 7 second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{7}))
xlim(xx)
ylim(yy)
subplot(3,2,4);hold on;plot(t,R);title('Quad muslce Left');
%  legend(fields{2},fields{11});title(trialgrf)
 xlabel('Time (s)');ylabel(' Muscle force (N)')
% legend(trials)
xlim(xx)
ylim(yy)


%hip left
[F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),LSol(locLast1:locLast2),100)
%all 7 second
% [F R]=ZeroTo100(GRFs(ind).data.(fields{end})(locLast1:locLast2),GRFs(ind).data.(fields{7}))
xlim(xx)
ylim(yy)
subplot(3,2,6);hold on;plot(t,R);title('Sol muslce Left');
%  legend(fields{2},fields{11});title(trialgrf)
 xlabel('Time (s)');ylabel(' Muscle force (N)')
% legend(trials)
xlim(xx)
ylim(yy)
end





 
%% Simple cognitive variables Performance:
% lets load the SH data from experiemnt and trigger etc and then plot the performance 
load('C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\Publications\Biomecognition_PaperMethod\Results\SH')
clear a b
subplot(1,2,1)
k=1;
for i = 1:length(SH.data.TrialNo)
if ~isnan(SH.data.Performance(i))
    a(k)= SH.data.Performance(i);
    b(k)= SH.data.TrialType(i);
k=k+1;
end
end
c=categorical(string(b));
bar(c,a)
xlabel('Trials');ylabel('Performance (%)')
ylim([0 100])

clear a b
subplot(1,2,2)
k=1;
for i = 1:length(SH.data.TrialNo)
if ~isnan(SH.data.HitPerformance(i))
    a(k)= SH.data.HitPerformance(i);
    b(k)= SH.data.TrialType(i);
k=k+1;
end
end
c=categorical(string(b));
bar(c,a)
xlabel('Trials');ylabel('Hit Performance (%)')
ylim([0 100])

%% •Target search performance

%% •Trigger Hit Performance

%% Cognitive motor interferences (CMI)


%% MoS dual tasking cost
clc
clear
load('C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\trcResults\gaitUtils.mat')
% W1W2W2 vs. W1+T... i.e. 17-18 (W1+T01-W1), 33-34, 49,50 in 
% gaitUtils.mos_AP(:,17) % this is W1+T
k=1;
for i=[18 34 50 17 33 49 ]
    a(k)=min(gaitUtils.mos_AP(:,i));
    b(k)=(gaitUtils.trialNames(:,i));
    k=k+1;
end
subplot(1,2,1)
c=categorical(string(b));
bar(c,a)
xlabel('Trials');ylabel('MoS: dual Tasking vs. ST in AP direction');
% 
% Cost of dual tasking = (ST MoS - DT MoS)/ST MoS X 100
% 
subplot(1,2,2)
for i=1:3
costMoS(i) = (a(i)-a(i+3))/(a(i))*100;
end
d=categorical(string(b(1:3)));
bar(d,costMoS)
xlabel('Trials');ylabel('Cost of dual Tasking on MoS in AP direction (%)');

%% •GRF

%% •Joint Kinematics and kinetics


%% •Muscle contribution

%% •Target search performance


%% •Cost of dual tasking


