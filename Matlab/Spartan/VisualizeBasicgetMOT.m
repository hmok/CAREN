
%% notes for future
% I added LG1Gait27_Reserve_Actuators.xml and everything works perfectly
% I also did filtering of GRF in Mokka, this is also needed to be done in
% Matlab and then all should be good.
% This is Hossein perfroming prelim study to bring C3D to Opensim
% and then to do Sclaing, IK, ID and SO so there are many issues so far. I
% am solving it step by step
% now (Mar 3rd 2019) it seems that GRF is not being read properly why in ID and this SO?
% howver, Scaling, IK are quite fine. A bit of fine tuning that may be done
% by students after the leanred the details. Let's focuson GRF, having 
% I think  I know waht the problem is. There are lots of Nan in mot file
% so Let's only focus on W1+L+ACC.c3d and start with no filtering for
% now. do the steps from c3d export to ID and SO.! scaling is done so we do
% IK onward. model: DFScaledv2.osim and all the nonneccessary joints are
% locked

% I am woring on this test folder for now:
% C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\TestOpenSim

% today 11:58am it seems ID is wokring. I am not doing any filtering so let's see
% how SO works now.

%% source folder
clear 
clc
close all
% source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH-Filtered';
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
tt=50;
grf = fullfile(source, 'grfResults'); trialgrf = strcat(trials{tt},'.mot');
trc = fullfile(source, 'trcResults');trialtrc = strcat(trials{tt},'.trc');%'W304.trc'
SO = fullfile(source, 'SOResults'); trialSO = strcat(trials{tt},'_StaticOptimization_force.sto');
ID = fullfile(source, 'IDResults');trialID = strcat(trials{tt},'_ID.sto');
IK = fullfile(source, 'IKResults');trialIK = strcat(trials{tt},'_IK.mot');

%% run Scale, IK, ID, SO, GRF dont seem to be working in ID and SO
% I locked mtp, subtalar, and lumbar DOFs be mindful of it
% naming: S = 0, T = 1, P = 2, W1 = 


% import org.opensim.modeling.*
% run('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab\c3dExport.m');
% !scale -S ScaleHBM_Setup.xml
% !ik -S IK_HBM_Setupv3.xml
% !opensim-cmd run-tool ID_HBM_Setup.xml 
% !opensim-cmd run-tool SO_HBM_Setupv3.xml

%% GRFs
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
figure
 hold on;plot(Forces.(fields{end}),Forces.(fields{2}));plot(Forces.(fields{end}),Forces.(fields{11}))
 legend(fields{2},fields{11});title(trialgrf)
%% trc plot
figure
import org.opensim.modeling.*
cd(fullfile(source,'trcResults'))

trctimeSeriesTable = TRCFileAdapter.read(trialtrc);
trc = osimTableToStruct(trctimeSeriesTable);
fields = fieldnames(trc);
figure
for i = 1:length(fields)-1
hold on
    plot(trc.(fields{i}))
end
% Left Foot 12:16, average X direction
figure
hold on;plot(trc.(fields{5}));plot(trc.(fields{6}))
legend(fields{5},fields{6})

% check if there is nan in the trc files
fnames = findnanfields(trc);



%% IK reslts plot
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



%% SO reslts plot
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

