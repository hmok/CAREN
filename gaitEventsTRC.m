%% finding the TO and HS using TRC files
% addpath('C:\OpenSim 4.0\Code\Matlab\');
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
import org.opensim.modeling.*
[path] = uigetdir('Select the folder that contain file.trc','*.trc');

%% let's loop for TRCs
for k=28%:2%

trcp = dir(fullfile(path,'*.trc'));
trcpath = fullfile(path,trcp(k).name)
trctimeSeriesTable = TRCFileAdapter.read(trcpath);
trc = osimTableToStruct(trctimeSeriesTable);
fields = fieldnames(trc);

    
% t=markerStruct.time.'; 
t = trc.time;
markerStruct.LMT2 = trc.LMT2;
markerStruct.RMT2 = trc.RMT2;
markerStruct.LHEE = trc.LHEE;
markerStruct.RHEE = trc.RHEE;
delta=30; %Parameter used for avoid peaks located closed to each other
% finding the toes
x_TL=markerStruct.LMT2(:,2); %Left Toe
x_TL=x_TL.';
x_TR=markerStruct.RMT2(:,2); %Right Toe
x_TR=x_TR.';
% finding the heels
x_HL=markerStruct.LHEE(:,2); %Left Heel
x_HL=x_HL.';
x_HR=markerStruct.RHEE(:,2); %  Right Heel
x_HR=x_HR.';
%Toes Max and Min Displacement
[maxtab_TL, mintab_TL] = peakdet_TL(x_TL, delta, t);
[maxtab_TR, mintab_TR] = peakdet_TR(x_TR, delta, t);
%Heels Max and Min Displacement
[maxtab_HL, mintab_HL] = peakdet_HL(x_HL, delta, t);
[maxtab_HR, mintab_HR] = peakdet_HR(x_HR, delta, t);


%% Plots

%% Plot compilation all 
figure(7) %-Y AXIS: Toe Left and Heel left upward and downward and forces Left and Right
markerStruct.time = t;
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

hold on; plot(mintab_HL(:,1), mintab_HL(:,2), 'g*','DisplayName','IC');
hold on

xlabel('Time (s)');
legend 
grid on
title('Left gait cycle');

%% X- Axis Analasys

%% finding the Toe and Heels
% t=markerStruct.time;
t=t.';
delta=10; %Parameter used for avoid peaks located closed to each other
% finding the toes
y_TL=markerStruct.LMT2(:,1); %Left Toe
y_TL=y_TL.';
y_TR=markerStruct.RMT2(:,1); %Right Toe
y_TR=y_TR.';
% finding the heels
y_HL=markerStruct.LHEE(:,1); %Left Heel
y_HL=y_HL.';
y_HR=markerStruct.RHEE(:,1); %  Right Heel
y_HR=y_HR.';
%Toes Max and Min Displacement
[maxtab_TL_x, mintab_TL_x] = peakdet_TL_x(y_TL, delta, t);
[maxtab_TR_x, mintab_TR_x] = peakdet_TR_x(y_TR, delta, t);
%Heels Max and Min Displacement
[maxtab_HL_x, mintab_HL_x] = peakdet_HL_x(y_HL, delta, t);
[maxtab_HR_x, mintab_HR_x] = peakdet_HR_x(y_HR, delta, t);

figure(8) %-Y AXIS: Toe Left and Heel left upward and downward and forces Left and Right

plot(markerStruct.time,markerStruct.LMT2(:,1),'b', 'LineWidth',2,'DisplayName','Left Toe 2')
hold on
plot(markerStruct.time,markerStruct.LHEE(:,1),'k', 'LineWidth',2,'DisplayName','Left Heel')
hold on

plot(markerStruct.time,markerStruct.RMT2(:,1),'--','DisplayName','Right Toe 2')
hold on
plot(markerStruct.time,markerStruct.RHEE(:,1),'--','DisplayName','Right Heel')
hold on
ylabel('Displacement [mm] X-axis');
%Left leg
plot(mintab_TL_x(:,1), mintab_TL_x(:,2), 'r*','DisplayName','MinT');
hold on; plot(mintab_HL_x(:,1), mintab_HL_x(:,2), 'g*','DisplayName','MinH');
plot(maxtab_TL_x(:,1), maxtab_TL_x(:,2), 'b*','DisplayName','MaxT');
hold on; plot(maxtab_HL_x(:,1), maxtab_HL_x(:,2), 'k*','DisplayName','MaxH');
%right leg
plot(mintab_TR_x(:,1), mintab_TR_x(:,2), 'r*','DisplayName','MinT');
hold on; plot(mintab_HR_x(:,1), mintab_HR_x(:,2), 'g*','DisplayName','MinH');
plot(maxtab_TR_x(:,1), maxtab_TR_x(:,2), 'b*','DisplayName','MaxT');
hold on; plot(maxtab_HR_x(:,1), maxtab_HR_x(:,2), 'k*','DisplayName','MaxH');

hold on

xlabel('Time (s)');
legend 
grid on
title('Left gait cycle');
%% Z- Axis Analasys

%% finding the Toe and Heels
t=markerStruct.time.';
% t=t.';
delta=10; %Parameter used for avoid peaks located closed to each other
% finding the toes
y_TL=markerStruct.LMT2(:,3); %Left Toe
y_TL=y_TL.';
y_TR=markerStruct.RMT2(:,3); %Right Toe
y_TR=y_TR.';
% finding the heels
y_HL=markerStruct.LHEE(:,3); %Left Heel
y_HL=y_HL.';
y_HR=markerStruct.RHEE(:,3); %  Right Heel
y_HR=y_HR.';
%Toes Max and Min Displacement
[maxtab_TL_z, mintab_TL_z] = peakdet_TL_z(y_TL, delta, t);
[maxtab_TR_z, mintab_TR_z] = peakdet_TR_z(y_TR, delta, t);
%Heels Max and Min Displacement
[maxtab_HL_z, mintab_HL_z] = peakdet_HL_z(y_HL, delta, t);
[maxtab_HR_z, mintab_HR_z] = peakdet_HR_z(y_HR, delta, t);

figure(9) %-z AXIS: Toe Left and Heel left upward and downward and forces Left and Right

plot(markerStruct.time,markerStruct.LMT2(:,3),'b', 'LineWidth',2,'DisplayName','Left Toe 2')
hold on
plot(markerStruct.time,markerStruct.LHEE(:,3),'k', 'LineWidth',2,'DisplayName','Left Heel')
hold on

plot(markerStruct.time,markerStruct.RMT2(:,3),'--','DisplayName','Right Toe 2')
hold on
plot(markerStruct.time,markerStruct.RHEE(:,3),'--','DisplayName','Right Heel')
hold on
ylabel('Displacement [mm] z-axis');
%Left leg
plot(mintab_TL_z(:,1), mintab_TL_z(:,2), 'r*','DisplayName','MinT');
hold on; plot(mintab_HL_z(:,1), mintab_HL_z(:,2), 'g*','DisplayName','MinH');
plot(maxtab_TL_z(:,1), maxtab_TL_z(:,2), 'b*','DisplayName','MaxT');
hold on; plot(maxtab_HL_z(:,1), maxtab_HL_z(:,2), 'k*','DisplayName','MaxH');
%right leg
plot(mintab_TR_z(:,1), mintab_TR_z(:,2), 'r*','DisplayName','MinT');
hold on; plot(mintab_HR_z(:,1), mintab_HR_z(:,2), 'g*','DisplayName','MinH');
plot(maxtab_TR_z(:,1), maxtab_TR_z(:,2), 'b*','DisplayName','MaxT');
hold on; plot(maxtab_HR_z(:,1), maxtab_HR_z(:,2), 'k*','DisplayName','MaxH');
hold on

xlabel('Time (s)');
legend 
grid on
% title('Left gait cycle');


end