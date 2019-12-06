%% ORS idea
% Title: stepping strategies during perturbed walking with/out congnitive task

%% 1. just focus on SH subject so source the folder
clear 
clc
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');

source ='C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';

%% 2. run this beaut below for MoS and other relvenat mearue of gait:
[mos_AP mos_ML bos_AP bos_ML, com_vel XCOM COM] = MoS(source,-1)
close all

%% 3. run the gaitEvents for step length, etc
 
 load gaitEvents.mat % from trcPath
%  figure;title('TO Left side');vline(gaitEvents(7).data.HSLocL) % plot the TO for left side or
%  figure;title('TO right side');vline(gaitEvents(7).data.HSLocR) % plot the TO for right side
%% plot step length for differen speeds and four trials
% trials = {'','','',''}

 pp =mean(mos_AP);
 ppstd =std(mos_AP);
  ppML =mean(mos_ML)
  ppMLstd =mean(mos_ML)
% g={'W1+L+ACC.trc';'W1+L+DEC.trc';'W1+R+ACC.trc';'W1+R+DEC.trc';'W1+T+L+ACC.trc';'W1+T+L+DEC.trc';'W1+T+R+ACC.trc';'W1+T+R+DEC.trc'}
% c= categorical(g)

% {  18     34     49        7         23        38         11          27            42         14       30       45 }
g={'W101','W201','W301','W1+L+DEC','W2+L+DEC','W3+L+DEC','W1+T+L+DEC','W2+T+L+DEC','W3+T+L+DEC','W1+T01','W2+T01','W3+T01' };
%     ;'W1+L+DEC.trc';'W1+R+ACC.trc';'W1+R+DEC.trc';'W1+T+L+ACC.trc';'W1+T+L+DEC.trc';'W1+T+R+ACC.trc';'W1+T+R+DEC.trc'}
c= categorical(g);

% bar(c,pp(6:13)) % first speed
% bar(c,pp(22:29))%seecond
% bar(c,pp(37:44))%third
subplot(1,2,1)
p3d = [ppML([18 34 49]);ppML([7 23 38]);ppML([11 27 42]);ppML([14 30 45])]
p3dstd = [ppMLstd([18 34 49]);ppMLstd([7 23 38]);ppMLstd([11 27 42]);ppMLstd([14 30 45])]
% bar3(p3d)

bar3_std(p3d,p3dstd)
set(gca,'XTickLabel',['W1'; 'W2';'W3'])
set(gca,'YTickLabel',['Trial1'; 'Trial2'; 'Trial3'; 'Trial4'])
zlabel('MoS (m)')
title('MoS in ML')
zlim([0 .6])
subplot(1,2,2)
p3d = [pp([18 34 49]);pp([7 23 38]);pp([11 27 42]);pp([14 30 45])]
p3dstd = [ppstd([18 34 49]);ppstd([7 23 38]);ppstd([11 27 42]);ppstd([14 30 45])]
% bar3(p3d)

bar3_std(p3d,p3dstd)
set(gca,'XTickLabel',['W1'; 'W2';'W3'])
set(gca,'YTickLabel',['Trial1'; 'Trial2'; 'Trial3'; 'Trial4'])
zlabel('MoS (m)')
title('MoS in AP')
zlim([0 .6])

