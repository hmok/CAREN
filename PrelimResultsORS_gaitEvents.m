
% this is to get MoS, gait paramers for static, dynamic
% question: what does pertubation do in different speed on gait speed
% 
% 
% L R ACC DEC
% first get the normal no pertub and differen speeds

% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\Trials')
% load all the trails names
close all
clear
clc
load('C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\Trials\trials.mat')

% e.g. trials{1,2} = 'S+T01'

% 5:8 = W1, 21:24 = W2,37:40 = W3 in trials so we need to find them in
% gaitEvents or any other results as the names are inculded. 

% lets talk about SH and find gaitEvent from trcFolder

source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH';

load(fullfile(source,'trcResults','gaitEvents.mat'))

% as for gaitEvents we have 18:21 = W1, 34:36 = W2,49:52 = W3 in trials so we need to find them in

% average them and plot them for all MoS, 
% strideLength

%% get the bar3(Z) e.g.:  load count.dat;Z = count(1:10,:);figure;bar3(Z);title('Detached Style')
% lets compare all for now from walking, pertub,T etc in 3D form, their
% average
% aligning axis of 3D use the following: https://github.com/phymhan/matlab-axis-label-alignment
% [L+ACC L+T+ACC T01 W0 ]
figure
j=1;
for i=[6 10 14 18]
sp(j,1)= mean(gaitEvents(i).data.stepLengthL);
j=j+1;
end

j=1;
for i=[22 26 30 34]
sp(j,2)= mean(gaitEvents(i).data.stepLengthL);
j=j+1;
end
j=1;
for i=[37 41 45 49 ]
sp(j,3)= mean(gaitEvents(i).data.stepLengthL);
j=j+1;
end
bar3(sp)
xlabel('gait Speed 2, 4 and 6km/hr');ylabel('Type of Trial, see title');zlabel('Step Length (m)')
title('1: Left ACC, 2: Left Target ACC,  3: Left Target Only, 4: Normal walking')

% https://github.com/phymhan/matlab-axis-label-alignment
%% some playing with the data to get a sense

% average left and right step length different speed

sp1= gaitEvents(18).data.stepLengthL;
sp2= gaitEvents(36).data.stepLengthL;
sp3= gaitEvents(49).data.stepLengthL; 
all = [mean(sp1);mean(sp2);mean(sp3)];
bar(all); title('Speed Effect');ylabel('Step Length(m)');
% legend('W1','W2','W3')

figure
sp1= gaitEvents(18).data.stepLengthR;
sp2= gaitEvents(36).data.stepLengthR;
sp3= gaitEvents(49).data.stepLengthR; 
all = [mean(sp1);mean(sp2);mean(sp3)];
bar(all); title('Speed Effect');ylabel('Step Length(m)');
% legend('W1','W2','W3')

% w1LACC =6 w2LACC=22, w3LACC =37
figure
sp1= gaitEvents(6).data.stepLengthL;
sp2= gaitEvents(22).data.stepLengthL;
sp3= gaitEvents(37).data.stepLengthL; 
all = [mean(sp1);mean(sp2);mean(sp3)];
bar(all); title('Left Leg, Speed Effect w/Pertubation on L and ACC');ylabel('Step Length(m)');
% legend('W1','W2','W3')

% w1RDEC =9 w2RDEC=25, w3RDEC =40
figure
sp1= gaitEvents(9).data.stepLengthR;
sp2= gaitEvents(29).data.stepLengthR;
sp3= gaitEvents(40).data.stepLengthR; 
all = [mean(sp1);mean(sp2);mean(sp3)];
bar(all); title('Right Leg, Speed Effect w/Pertubation on R and DEC');ylabel('Step Length(m)');
% legend('W1','W2','W3')