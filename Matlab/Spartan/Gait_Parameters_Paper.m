%% Paper 1, Gait parameters due to pertubation L/R, Acc/Dec
% this is for the paper Gait parameters and 
% There are several ideas I can write from the below:

% 	Parameters: gait speeds×target search ×perturbation
% 	Parameters: gait speeds×perturbation [current paper]
% 	Parameters: target search ×perturbation
% 	Parameters: gait speeds×target search

%% the focus in this study is the follwoing i.e. 15:
% W1-3 and W1+P - W3+P , P=L/R and ACC/DEC
trials = {'W1', 'W1+L+ACC',  'W1+L+DEC',  'W1+R+ACC',  'W1+R+DEC', ...
    'W2', 'W2+L+ACC',  'W2+L+DEC', 'W2+R+ACC',  'W2+R+DEC', ...
    'W3', 'W3+L+ACC',  'W3+L+DEC', 'W3+R+ACC',  'W3+R+DEC'};

%% for simplicity just focus on W1-3 and L+ACC i.e. 6 trials

% Though I am focuding on the second idea i.e.	Parameters: gait speeds×perturbation [current paper]
% Just focus on one person and get all the results and then add more people
% So now I am focusing on SH in this folder: source + SH
% source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH';

%% Author(s): Hossein Mokhtarzadeh
% v1 26 Aug 2019
% To be submitted in 30  Oct 2019 to Journal of Biomechanics

%% Location
% C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\Publications\Gait_Parameters_Paper


%% Figure 1
% Explain: peak GRF (even rate to peak as well maybe) before and after pertubation
% average standing vs. Normal vs. perturbed using bar
% so far I cannto undersntad what this can help but not bad to get an idea of max GRF in diffrente senariios
clear
close all
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab')
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH';
cd(fullfile(source,'grfResults'))
load('GRFs.mat')
k=1;
% for i=[2,6,7,8,9,18,22:25,34,37,38,40,49 ]%2:length(GRFs)%[2,6,18,22,34,37,49]
for i=[6,8,18,22, 24,34,37,39,49 ]%2:length(GRFs)%[2,6,18,22,34,37,49]
maxes(k)  = max(GRFs(i).data.ground_force_1_vy);
c{k} = GRFs(i).name;
k=k+1;
end
% d= categorical(c);
figure
bar(maxes)
set(gca, 'XTickLabel', c)
xtickangle(45)
title('Max GRF in 7s')
ylabel('Peak GRF (N)')

%% Figure 1a, I dpont like these graphs as there are errors in CoP
% Explain: peak GRF COP (even rate to peak as well maybe) before and after pertubation
% average standing vs. Normal vs. perturbed using bar
% so far I cannto undersntad what this can help but not bad to get an idea of max GRF in diffrente senariios
clear
close all
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH';
cd(fullfile(source,'trcResults'))
load('gaitEvents.mat')
cd(fullfile(source,'grfResults'))
load('GRFs.mat')

k=1;
for i=[6,8,18,22, 24,34,37,39,49 ]%2:length(GRFs)%[2,6,18,22,34,37,49]
maxes(k)  = max(GRFs(i).data.ground_force_1_px); % be casreful of large CoP
c{k} = GRFs(i).name;
k=k+1;
end
% d= categorical(c);
figure
bar(maxes)
set(gca, 'XTickLabel', c)
xtickangle(45)
title('Max CoP in 7s')
ylabel('Peak CoP (m)')

figure
trial = 28;
[I J1]=min(abs(gaitEvents(trial).data.HSLocL(1)-GRFs(trial).data.time))
[I J2]=min(abs(gaitEvents(trial).data.HSLocL(2)-GRFs(trial).data.time))
x=(GRFs(trial).data.ground_force_1_px(J1:J2)); % be casreful of large CoP
y=(GRFs(trial).data.ground_force_1_pz(J1:J2)); % be casreful of large CoP

plot(x,y);xlabel('xCoP');ylabel('zCoP');title(GRFs(trial).name)
gaitEvents(trial).name

% xlim([-1 1]);ylim([-1 1])

%% Figure 2
% Explain: MoS and other gait parameters bos, CoM vel, CoM, etc
% Let's do mos in AP, ML
clear
close all
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH';
cd(fullfile(source,'trcResults'))
load('gaitUtils.mat')

% for i=[2,4]%[2,6,7,8,9,18,22:25,34,37,38,40,49 ]%2:length(GRFs)%[2,6,18,22,34,37,49]
minsAP  = mean(gaitUtils.mos_AP);%min(gaitUtils.mos_AP);
minsML  = mean(gaitUtils.mos_ML);%min(gaitUtils.mos_ML);
c = gaitUtils.trialNames;

% d= categorical(c);
show=[6:9,18,22:25,34,37:40,49]
figure
bar(minsAP(:,show))
set(gca, 'XTickLabel',c(show))
xtickangle(45)
title('MoS AP')

figure
bar(minsML(:,show))
set(gca, 'XTickLabel',c(show))
xtickangle(45)
title('MoS ML')
%% Figure 3
% Explain: Get gaitevents step length and step freq etc
clear
close all
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH';
cd(fullfile(source,'trcResults'))
load('gaitEvents.mat')

k=1;
show=[6:9,18,22:25,34,37:40,49];
for i=show%[2,6,7,8,9,18,22:25,34,37,38,40,49 ]%2:length(GRFs)%[2,6,18,22,34,37,49]
meansSL(:,k)  = mean(gaitEvents(i).data.stepLengthL);
% minsML  = min(gaitUtils.mos_ML);
c{k} = gaitEvents(i).name;
k=k+1;
end
% d= categorical(c);

figure
bar(meansSL)
set(gca, 'XTickLabel',c)
xtickangle(45)
title('Step Length (Left)')
ylabel('Step Length (m)')

%% Figure 4
% Explain: Get gaitevents step width Left side
clear
close all
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH';
cd(fullfile(source,'trcResults'))
load('gaitEvents.mat')

k=1;
show=[6:9,18,22:25,34,37:40,49];
for i=show%[2,6,7,8,9,18,22:25,34,37,38,40,49 ]%2:length(GRFs)%[2,6,18,22,34,37,49]
meansSWL(:,k)  = mean(gaitEvents(i).data.stepWidthL);
% minsML  = min(gaitUtils.mos_ML);
c{k} = gaitEvents(i).name;
k=k+1;
end
% d= categorical(c);

figure
bar(meansSWL)
set(gca, 'XTickLabel',c)
xtickangle(45)
title('Step Width (Left)')
ylabel('Step Width (m)')

%% Figure 5
% Explain: Get gaitevents step width Right side
clear
close all
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH';
cd(fullfile(source,'trcResults'))
load('gaitEvents.mat')

k=1;
show=[6:9,18,22:25,34,37:40,49];
for i=show%[2,6,7,8,9,18,22:25,34,37,38,40,49 ]%2:length(GRFs)%[2,6,18,22,34,37,49]
meansSWR(:,k)  = mean(gaitEvents(i).data.stepWidthR);
% minsML  = min(gaitUtils.mos_ML);
c{k} = gaitEvents(i).name;
k=k+1;
end
% d= categorical(c);

figure
bar(meansSWR)
set(gca, 'XTickLabel',c)
xtickangle(45)
title('Step Width (Right)')
ylabel('Step Width (m)')


%% Figure 6
% Explain: Get gaitevents step frqucny Left side
clear
close all
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH';
cd(fullfile(source,'trcResults'))
load('gaitEvents.mat')

k=1;
show=[6:9,18,22:25,34,37:40,49];
for i=show%[2,6,7,8,9,18,22:25,34,37,38,40,49 ]%2:length(GRFs)%[2,6,18,22,34,37,49]
meansSFL(:,k)  = mean(gaitEvents(i).data.strideFreqL);
% minsML  = min(gaitUtils.mos_ML);
c{k} = gaitEvents(i).name;
k=k+1;
end
% d= categorical(c);

figure
bar(meansSFL)
set(gca, 'XTickLabel',c)
xtickangle(45)
title('Stride Frequency (Left)')
ylabel('Stride Frequency (1/s)')


%% Figure 7: Perturbation effect 
% Explain: So far we had the average effect how about the role of
% peturbation on these values e.g. MoS, SL, etc before and after perturb 
% Let's jsut focus on W1 L+ACC pre-post pertubr

clear
close all
trials = {'W101', 'W1+L+ACC',  'W1+L+DEC',  'W1+R+ACC',  'W1+R+DEC', ...
    'W201', 'W2+L+ACC',  'W2+L+DEC', 'W2+R+ACC',  'W2+R+DEC', ...
    'W301', 'W3+L+ACC',  'W3+L+DEC', 'W3+R+ACC',  'W3+R+DEC'};

addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab')
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH';
cd(fullfile(source,'DFLOW'))% LACC = 28
load('dflow.mat')
iniP=dflow.Perturb(26).iniPerturb;endP=dflow.Perturb(26).endPertub;
cd(fullfile(source,'trcResults'))
load('gaitUtils.mat')
for i=1:length(dflow.Perturb)
 trials2{i}= dflow.Perturb(i).name
end

% Match the index of trials and names from dflow

[ind missTrials] = trialIndexFinder(trials,trials2);

% for i=[2,4]%[2,6,7,8,9,18,22:25,34,37,38,40,49 ]%2:length(GRFs)%[2,6,18,22,34,37,49]
mosallAP  = mean(gaitUtils.mos_AP);%min(gaitUtils.mos_AP);
mosallML  = mean(gaitUtils.mos_ML);%min(gaitUtils.mos_ML);
mosPreAP  = mean(gaitUtils.mos_AP(1:iniP,:));%min(gaitUtils.mos_AP);
mosPreML  = mean(gaitUtils.mos_ML(1:iniP,:));%min(gaitUtils.mos_ML);
mosPosAP  = mean(gaitUtils.mos_AP(endP:end,:));%min(gaitUtils.mos_AP);
mosPosML  = mean(gaitUtils.mos_ML(endP:end,:));%min(gaitUtils.mos_ML);

c = gaitUtils.trialNames;

% d= categorical(c);
show=[6 ,18];
figure
bar([mosallAP(:,6); mosallAP(:,18); mosPreAP(:,6) ;mosPosAP(:,6)])
d={'mosallAPW1LACC','mosallAPW1','mosPreAPW1LACC','mosPosAPW1LACC'};
d={'All-Pert','All-W1','Pre-Pert','Pos-Pert'};
% bar(mosallAP(:,show))
% bar(minsAP(:,show))
% set(gca, 'XTickLabel',c(show))
set(gca, 'XTickLabel',d)
xtickangle(45)
title('MoS AP')

 figure
bar([mosallML(:,6); mosallML(:,18); mosPreML(:,6) ;mosPosML(:,6)])
 % bar(minsML(:,show))
% set(gca, 'XTickLabel',c(show))
 set(gca, 'XTickLabel',d)
 xtickangle(45)
 title('MoS ML')


%% Table 1
% Explain: create a table for all 15 trails this way, their averages and
% sd
% mos,bos from gaitUtils, stepL,Freq,Width from gaitEvents

%             | MoSAP | MoSML | StepL | StepFreq | bosAP | bosML | vGRF |
% |    W1    | mean(sd)
% |    W2    |
% |    W3    |
% |    ..    |    
% | W3+R+DEC |
% Let's find the means for now

% Average MoS in AP for W1

% Average StepLength for W1 (average L+R)



%% Table 2
% Explain:

%% Table 3
% Explain:

%% Table 4
% Explain: