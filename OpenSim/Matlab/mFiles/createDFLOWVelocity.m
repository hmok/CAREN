%% to create the following velocity profile using parameters we input in CAREN
%first version I made in July 5th 2019
% TODOs: create a proper function 
% add left/right and make it more generalized.
% check nargin == 6 as we need 6 points to recreate the profile.

%% function [v] = createDFLOWVelocity(offset,Acc,Time, w1,t,tP,tA,freq,ploting)

% States & Variables:

%              |  Time  |
%     -        ..........
%             .          .
%   offset   .            . 
%           .              .
% ..........                ...........
% | 0  | 1 | 2 |   3    | 4 |   0       These are the different states per DoF / Belt


% Slope is determined by Acc

offset = .6;
Acc = 0.3;
Time = .2;
w1 = 1.667;
t = 7; %t = tP+tA;
tP = 2; %time before pertubr
tA = 5;%time afer pertubr
interval = 0.0033;%freq = 300;

p1=[w1,0]; %the points 6 of them to create DFLOW velocity 
p2=[w1,2]; %the points 6 of them to create DFLOW velocity 
p3=[w1+offset,Acc/w1]; %the points 6 of them to create DFLOW velocity 
p4=[w1+offset,Acc/w1+Time]; %the points 6 of them to create DFLOW velocity 
p5=[w1,2*Acc/w1+Time]; %the points 6 of them to create DFLOW velocity 
p6=[w1,7]; %the points 6 of them to create DFLOW velocity 

y =[w1;w1;w1+offset;w1+offset;w1;w1];
x =[0;2;2+Acc/w1;2+Acc/w1+Time;2+2*Acc/w1+Time;7];
time = 0:interval:t;
plot(x,y,'o',time,[pchip(x,y,time)])
legend('data','pchip','Location', 'Best')
%     legend('data','pchip','spline', 'Location', 'Best')

v = pchip(x,y,time);
plot(time,v)
[Xa,Ya] = alignsignals(y,v)
% plot(x,y,'o',time,[pchip(x,y,time); spline(x,y,time)])
% 
% a=Data(:,RbeltSpeed)
% plot(a)
% find max and min
% % findpeaks(a,'MinPeakProminence',40)
% 
% %  findpeaks(a,'MinPeakHeight',0.1,'MinPeakDistance',20)
% [maxtab, mintab]=peakdet_TL(a, delta, time)
% delta = 30;
% peakdet(a,delta,time)
% peakdet_TL(a, delta, time)
% plot(time,a)
