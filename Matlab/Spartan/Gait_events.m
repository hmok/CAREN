clc
close all
% clear all
% Import *.c3d file inside Matlab
% addpath('C:\OpenSim 4.0\Code\Matlab\');
% c3dExport


%% Analysis 
%% finding the min and max peaks
t=markerStruct.time; % import(trc)
t=t.';
delta=45;
% findin the toes
x_TL=markerStruct.LMT2(:,2);
x_TL=x_TL.';
x_TR=markerStruct.RMT2(:,2);
x_TR=x_TR.';
% findin the heels
x_HL=markerStruct.LHEE(:,2);
x_HL=x_HL.';
x_HR=markerStruct.RHEE(:,2);
x_HR=x_HR.';

[maxtab_TL, mintab_TL] = peakdet_TL(x_TL, delta, t);
[maxtab_TR, mintab_TR] = peakdet_TR(x_TR, delta, t);
%Heels
[maxtab_HL, mintab_HL] = peakdet_HL(x_HL, delta, t);
[maxtab_HR, mintab_HR] = peakdet_HR(x_HR, delta, t);
%% Plots

figure %-Y AXIS: Toe Left and Heel left upward and downward and forces Left and Right

for i=1:1:length(mintab_HL) %or length(mintab_HL) % positon must be greater>0 so i use abs but talk to Carlos
% rectangle('Position',abs([mintab_HL(i,1) 0 (mintab_TL(i,1)-mintab_HL(i,1)) mintab_HL(i,2)]),'FaceColor','m') %Stance phase
 hold on
% hLHS = vline(mintab_HL(i,1),'-r','LHS');
% hRHS = vline(mintab_HR(i,1),'r','RHS');
% hLTO = vline(mintab_TL(i,1),'-b','LTO');
% hRTO = vline(mintab_TR(i,1),'b','RTO');
 % 
% rectangle('Position',abs([mintab_TL(i,1) 0 (mintab_HL(i+1,1)-mintab_TL(i,1)) mintab_HL(i,2)]),'FaceColor','c') %Swing phase
hold on
end


plot(0,'Marker','square','LineStyle','none','Color','b','DisplayName','Stance phase','MarkerFaceColor','m');
hold on
plot(0,'Marker','square','LineStyle','none','Color','b','DisplayName','Swing phase','MarkerFaceColor','c');
hold on

plot(markerStruct.time,markerStruct.LMT2(:,2),'b', 'LineWidth',2,'DisplayName','Left Toe 2')
hold on
plot(markerStruct.time,markerStruct.LHEE(:,2),'k', 'LineWidth',2,'DisplayName','Left Heel')
hold on

plot(markerStruct.time,markerStruct.RMT2(:,2),'--','DisplayName','Right Toe 2')
hold on
plot(markerStruct.time,markerStruct.RHEE(:,2),'--','DisplayName','Right Heel')
hold on
ylabel('Displacement [mm]');

plot(mintab_TL(:,1), mintab_TL(:,2), 'r*','DisplayName','TO');
% plot(maxtab_TL(:,1), maxtab_TL(:,2), 'r*');
% hold on; plot(mintab_TR(:,1), mintab_TR(:,2), 'g*');
% plot(maxtab_TR(:,1), maxtab_TR(:,2), 'r*');

hold on; plot(mintab_HL(:,1), mintab_HL(:,2), 'g*','DisplayName','IC');
% plot(maxtab_HL(:,1), maxtab_HL(:,2), 'r*');
% hold on; plot(mintab_HR(:,1), mintab_HR(:,2), 'g*');
% plot(maxtab_HR(:,1), maxtab_HR(:,2), 'r*');

% end of finding peaks 
%%
% Force Right ad Left
yyaxis right
plot(forceStruct.time,forceStruct.f1(:,2),'k','DisplayName','G.Force Left')
hold on
plot(forceStruct.time,forceStruct.f2(:,2),'k','DisplayName','G.Force Right')
grid on
ylabel('Force Y-Axis [N]');

xlabel('Time (s)');
legend 
grid on
title('Left gait cycle');

IC = round(mintab_HL(:,1)*100+1);
TO = round(mintab_TL(:,1)*100+1);

%% IC and TO for both feet
for i = 1:length(IC)
hLHS = vline(IC(i),'-r','LHS');
hRHS = vline(mintab_HR(i,1),'r','RHS');
end
for i = length(IC)
hLHS = vline(,'-r','LHS');
hRHS = vline(mintab_HR(i,1),'r','RHS');
end
hLTO = vline(mintab_TL(i,1),'-b','LTO');
hRTO = vline(mintab_TR(i,1),'b','RTO');


%% Stance phase
for i=1:1:length(IC)
stance_LT{i} = IC(i,1)-1:1:TO(i,1)-1;
stance_LH{i} = IC(i,1):1:TO(i,1);
end

figure

for i=1:1:length(IC)

subplot (1,length(IC),i)
title(sprintf('Stance phase %i',i));
hold on
plot(stance_LT{i}/100,markerStruct.LMT2(stance_LT{i},2),'b', 'LineWidth',2,'DisplayName','Left Toe 2')
hold on
plot(stance_LH{i}/100,markerStruct.LHEE(stance_LH{i},2),'r', 'LineWidth',2,'DisplayName','Left Heel')
hold on
yyaxis right
plot(stance_LH{i}/100,forceStruct.f1(stance_LH{i}*20,2),'k','DisplayName','G.Force Left')
hold on
plot(stance_LH{i}/100,forceStruct.f2(stance_LH{i}*20,2),'k','DisplayName','G.Force Right')
hold on
grid on

end
legend 


%% swing phase
for i=1:1:length(IC)-1
swing_LT{i} = TO(i,1):1:IC(i+1,1);
swing_LH{i} = TO(i,1):1:IC(i+1,1);
end

figure

for i=1:1:length(IC)-1
subplot (1,length(IC)-1,i)

plot(swing_LT{i}/100,markerStruct.LMT2(swing_LT{i},2),'b', 'LineWidth',2,'DisplayName','Left Toe 2')
hold on
plot(swing_LH{i}/100,markerStruct.LHEE(swing_LH{i},2),'r', 'LineWidth',2,'DisplayName','Left Heel')
hold on
yyaxis right
plot(swing_LT{i}/100,forceStruct.f1(swing_LT{i}*20,2),'k','DisplayName','G.Force Left')
hold on
plot(swing_LT{i}/100,forceStruct.f2(swing_LT{i}*20,2),'k','DisplayName','G.Force Right')
hold on
grid on
title(sprintf('Swing phase %i',i));
end
legend 
