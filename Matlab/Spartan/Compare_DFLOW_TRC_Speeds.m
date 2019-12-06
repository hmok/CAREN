%% let's compare DFLOW and TRC belt speeds

% 0. get the folder needed...path and session , i may need to which trial in this session the list of trials are needed..
% 1. recreate belt speed and save using dflowReproduceTRC.m: change the dirctory or make it as function toget path
% 2. load BeltVel for this person sesseion e.g.,
% load('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\DFLOW\BeltVel.mat')
% {'W1+L+ACC.trc'}; figure;plot(BeltVel(8).data);hold on;plot(Lvel(:,end))
% 3. load actual belt speed if exist from the same session e,g, using dflow.m it asks for the folder the txt files are:
% so run this "dflow.m" dflow by the way, plots the relevant traisl
% dflow can provide all four traisl velocities for L/R i.e. Lvel(first:fourth) and Rvel(first:fourth)
% 
% 4. plot them how?

clc
clear
%% run them as below:
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab')
dflow
cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\trcResults')
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB2\VICON\CB2\trcResults')
% load TRCs
load gaitEvents.mat
load BeltVel.mat; %from dflowReproduceTRC_Events or dflowReproduceTRC
hold on; plot(BeltVel(29).data);title(BeltVel(29).name)
vline(gaitEvents(7).data.HSLocL*100);

% figure
% vline(gaitEvents(14).data.HSLocL*100);
% hold on; plot(BeltVel(14).data);
% hold on; plot(BeltVel(8).data);
% hold on; plot(BeltVel(9).data);
% hold on; plot(BeltVel(10).data);
% hline(0)
