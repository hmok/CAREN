% ************************************************************
% Create Structure Array
% https://www.mathworks.com/help/matlab/matlab_prog/create-a-structure-array.html
% Author(s): Hossein Mokhtarzadeh
% Date: 27 June 2019

% Describtion: Let's get all the TRC, Mot, MoS, etc files and get one Structure file i.e. TRCs.mat
% for visualiaton later, plot,etc
% ************************************************************


% TODOs: Let's visualize what I collected from Data_Collection_OpenSim.m
% Ask a question and then bring the geberal forms here to visualize
% but for a specific paper create a relevant M file so that I can refer to
% later


%% For visuliazaton later on
% Question: Lets compare GRF for DF1 (W1+T+P) i.e. [10:13] or GRFs(10:13)
% {'W1+T+L+ACC.mot';'W1+T+L+DEC.mot';'W1+T+R+ACC.mot';'W1+T+R+DEC.mot'}
% fields = fieldnames(GRFs(1).data)
% plot(GRFs(1).data.(fields{1}))
% GRFs(1).data
%let
load GRFs %load GRFs from the grfResults of its folder DF1
ToPlot = GRFs(10:13); %select the data
fields = fieldnames(ToPlot(1).data); %exmaple how to get fields
figure(1)
for i=1:length(ToPlot) %
    hold on;
    h(i) = plot(ToPlot(i).data.(fields{end}),ToPlot(i).data.(fields{2}));% only the left side
    names{i} = ToPlot(i).name;
end
legend(h, names)
% xlim([0 2])
figure(2)
for i=1:length(ToPlot) %
    hold on;
    h(i) = plot(ToPlot(i).data.(fields{end}),ToPlot(i).data.(fields{11}));% only the right side
    names{i} = ToPlot(i).name;
end
legend(h, names)

%% IK visualize
clear;clc
cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\IKResults');
load IKs %load GRFs from the grfResults of its folder DF1
ToPlot = IKs(3:end); %select the data
fields = fieldnames(ToPlot(1).data); %exmaple how to get fields
figure(1)
for i=1:length(ToPlot) %
    hold on;
    %do this later for accuate nameing etc: [a b]=ismember('ground_force_1_vy',fields) then b =2 is fields{2}
    h(i) = plot(ToPlot(i).data.(fields{end}),ToPlot(i).data.(fields{2}));% only the left side or do this later
    names{i} = ToPlot(i).name;
end
legend(h, names)
% xlim([0 2])
clear h; close all

% IK outcomes 
firstPlots = {'W1+L+ACC_ik.mot'}%,'W2+L+ACC_ik.mot','W2+L+ACC_ik.mot'};
secondPlots = {'W1+R+ACC_ik.mot'}%,'W2+R+ACC.mot_ik','W2+R+ACC_ik.mot'};
thirdPlots = {'W101_ik.mot'}%,'W2+L+ACC_ik.mot','W3+L+ACC_ik.mot'};
fourthPlots = {'W1+R+DEC_ik.mot'}%,'W2+R+DEC_ik.mot','W3+R+DEC_ik.mot'};

subplot(2,2,1) % left - ACC
title('Left-ACC')
for i=1:length(firstPlots) %
    hold on;
    [lia,locb] = ismember(firstPlots{i},names);
%     ind = find(a, 1, 'first');
ind = locb;
    h(i) = plot(ToPlot(ind).data.(fields{end}),ToPlot(ind).data.(fields{2}));% only the left side
%     names{i} = ToPlot(i).name;
    
end
legend(h, firstPlots)

subplot(2,2,2) % Right - ACC
title('Right-ACC')
for i=1:length(secondPlots) %
    hold on;
    [lia,locb] = ismember(secondPlots{i},names);
%     ind = find(a, 1, 'first');
ind = locb;
    h(i) = plot(ToPlot(ind).data.(fields{end}),ToPlot(ind).data.(fields{2}));% only the left side
%     names{i} = ToPlot(i).name;
    
end
legend(h, secondPlots)

subplot(2,2,3) % left - DEC
title('Left-DEC')
for i=1:length(thirdPlots) %
    hold on;
    [lia,locb] = ismember(thirdPlots{i},names);
%     ind = find(a, 1, 'first');
ind = locb;
    h(i) = plot(ToPlot(ind).data.(fields{end}),ToPlot(ind).data.(fields{2}));% only the left side
%     names{i} = ToPlot(i).name;
    
end
legend(h, thirdPlots)

subplot(2,2,4) % Right - DEC
title('Right-DEC')
for i=1:length(fourthPlots) %
    hold on;
    [lia,locb] = ismember(fourthPlots{i},names);
%     ind = find(a, 1, 'first');
ind = locb;
    h(i) = plot(ToPlot(ind).data.(fields{end}),ToPlot(ind).data.(fields{2}));% only the left side
%     names{i} = ToPlot(i).name;
    
end
legend(h, fourthPlots)


% figure(2)
% for i=1:length(ToPlot) %
%     hold on;
%     h(i) = plot(ToPlot(i).data.(fields{end}),ToPlot(i).data.(fields{11}));% only the right side
%     names{i} = ToPlot(i).name;
% end
% legend(h, names)
%%

clear;clc
cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\grfResults');
load GRFs %load GRFs from the grfResults of its folder DF1
ToPlot = GRFs(3:end); %select the data
fields = fieldnames(ToPlot(1).data); %exmaple how to get fields
figure(1)
for i=1:length(ToPlot) %
    hold on;
    %do this later for accuate nameing etc: [a b]=ismember('ground_force_1_vy',fields) then b =2 is fields{2}
    h(i) = plot(ToPlot(i).data.(fields{end}),ToPlot(i).data.(fields{2}));% only the left side or do this later
    names{i} = ToPlot(i).name;
end
legend(h, names)
% xlim([0 2])
clear h; close all
% first Left foot and Acceleation (W1 W2 W3) 2x2 table
firstPlots = {'S+T01.mot','W101.mot','W201.mot','W301.mot','W1+L+ACC.mot','W2+L+ACC.mot','W3+L+ACC.mot'};
% firstPlots = {'W1+L+ACC.mot','W2+L+ACC.mot','W3+L+ACC.mot'};
secondPlots = {'S+T01.mot','W101.mot','W201.mot','W301.mot','W1+L+DEC.mot','W2+L+DEC.mot','W3+L+DEC.mot'};
thirdPlots = {'S+T01.mot','W101.mot','W201.mot','W301.mot','W1+R+ACC.mot','W2+R+ACC.mot','W3+R+ACC.mot'};
fourthPlots = {'S+T01.mot','W101.mot','W201.mot','W301.mot','W1+R+DEC.mot','W2+R+DEC.mot','W3+R+DEC.mot'};

firstPlots = {'W1+L+ACC.mot','W2+L+ACC.mot','W3+L+ACC.mot'};
secondPlots = {'W1+L+DEC.mot','W2+L+DEC.mot','W3+L+DEC.mot'};
thirdPlots = {'W1+R+ACC.mot','W2+R+ACC.mot','W3+R+ACC.mot'};
fourthPlots = {'W1+R+DEC.mot','W2+R+DEC.mot','W3+R+DEC.mot'};

firstPlots = {'W1+L+ACC.mot','W2+L+ACC.mot','W2+L+ACC.mot'};
secondPlots = {'W1+R+ACC.mot','W2+R+ACC.mot','W2+R+ACC.mot'};
thirdPlots = {'W1+L+ACC.mot','W2+L+ACC.mot','W3+L+ACC.mot'};
fourthPlots = {'W1+R+DEC.mot','W2+R+DEC.mot','W3+R+DEC.mot'};

% IK outcomes 
% firstPlots = {'W1+L+ACC_ik.mot'}%,'W2+L+ACC_ik.mot','W2+L+ACC_ik.mot'};
% secondPlots = {'W1+R+ACC_ik.mot'}%,'W2+R+ACC.mot_ik','W2+R+ACC_ik.mot'};
% thirdPlots = {'W101_ik.mot'}%,'W2+L+ACC_ik.mot','W3+L+ACC_ik.mot'};
% fourthPlots = {'W1+R+DEC_ik.mot'}%,'W2+R+DEC_ik.mot','W3+R+DEC_ik.mot'};

subplot(2,2,1) % left - ACC
title('Left-ACC')
for i=1:length(firstPlots) %
    hold on;
    [lia,locb] = ismember(firstPlots{i},names);
%     ind = find(a, 1, 'first');
ind = locb;
    h(i) = plot(ToPlot(ind).data.(fields{end}),ToPlot(ind).data.(fields{2}));% only the left side
%     names{i} = ToPlot(i).name;
    
end
legend(h, firstPlots)

