%% grfGaitEvent(grfpath)
% function [LTO, LFS, RTO, RFS, timeLTO, timeLFS,timeRTO, timeRFS] = grfGaitEvent(RFy, LFy,threshold, ploting)
function gaitEvents = grfGaitEvent(RFy, LFy, t, threshold, ploting)
%anotehr idea is to just get GRFs and then run all [LTO, LFS, RTO, RFS] = grfGaitEvent(RFy, LFy,trial,threshold, ploting)

%% example:
% what does it do? this code provides all the timing,  gait events (frames of GRFs) given left or right vertical GRFs in CAREN
% How to use it? This needs four input (Right and Left forces, threshold=20 for cutoff, whether to plot or not=1)

%%%%%%% For instance run the following lines to get the gait events %%%%%%

%%%%%%% Begin %%%%%%
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\grfResults');
% load GRFs
% addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab')
% ToPlot = GRFs(3:end); %select the data
% fields = fieldnames(ToPlot(1).data); %exmaple how to get fields
% tr=25; %wchic trial
% RFy = ToPlot(tr).data.(fields{2}); 
% LFy = ToPlot(tr).data.(fields{11});
% % [LTO LFS RTO RFS] = grfGaitEvent(RFy, LFy,20, 1)
% time = ToPlot(tr).data.(fields{end}); %or t in the function
% gaitEvents = grfGaitEvent(RFy, LFy,time, 20, 1)
% time = ToPlot(tr).data.(fields{end}); %or t in the function
% timeLTO  = time(LTO);
% timeLFS  = time(LFS);
% timeRTO  = time(RTO);
% timeRFS  = time(RFS);
%%%%%%% End %%%%%%

% [LTO LFS RTO RFS] = grfGaitEvent(RFy, LFy,threshold=20, ploting=1)
% TODOs
% Pertubr is still challengin
% I think we need to identify pertub timing and then give enough time that
% delta (FS-FS) gets back to normal...
% how to use for large number of files
%% Main part of script

for i=1:4
    disp('be careful about the perturbed trials')
end
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\grfResults');
% load GRFs
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\trcResults');
% load TRCs


time = 1:length(RFy);% ToPlot(i).data.(fields{end});

y2 = threshold*ones(length(RFy),1);
[X0,Y0,I,J] = intersections(time,RFy,time,y2);

subplot(1,2,2)
plot(time, RFy);hold on;vline(X0)

% Cutoff = 20;
X0 = round(X0);
yLimits = get(gca,'YLim');
% right foot off and strike

k=1;j=1;
for i=1:numel(X0)
    if RFy(X0(i)+40) > RFy(X0(i)) % 40 more frames are greater than the current (GRF is increasing)
        RFS(j,1) = X0(i);
        text(RFS(j,1),yLimits(2)+30, {'RFS'},'VerticalAlignment', 'top')
        j=j+1;
    else
        RTO(k,1) = X0(i);
        text(RTO(k,1),yLimits(2)+30, {'RTO'},'VerticalAlignment', 'top')
        k=k+1;
    end
    
end

% left foot off and strike

% LFy = ToPlot(i).data.(fields{11});
time = 1:length(LFy);% ToPlot(i).data.(fields{end});

y2 = threshold*ones(length(LFy),1);
[X0,Y0,I,J] = intersections(time,LFy,time,y2);
% figure;title('Left foot')
subplot(1,2,1)
plot(time, LFy);hold on;vline(X0)

% Cutoff = 20;
X0 = round(X0);
yLimits = get(gca,'YLim');


% left foot off and strike
k=1;j=1;
for i=1:numel(X0)
    if LFy(X0(i)+40) > LFy(X0(i))
        LFS(j,1) = X0(i);
        text(LFS(j,1),yLimits(2)+30, {'LFS'},'VerticalAlignment', 'top')
        j=j+1;
    else
        LTO(k,1) = X0(i);
        text(LTO(k,1),yLimits(2)+30, {'LTO'},'VerticalAlignment', 'top')
        k=k+1;
    end
end


if ploting~=1
    close all
end
% SaveName = 'LR_HS'

timeLTO  = t(LTO);
timeLFS  = t(LFS);
timeRTO  = t(RTO);
timeRFS  = t(RFS);
I;J;Y0;
events = {LTO, LFS, RTO, RFS, timeLTO, timeLFS,timeRTO, timeRFS};
eventNames = {'LTO', 'LFS', 'RTO', 'RFS', 'timeLTO', 'timeLFS','timeRTO', 'timeRFS'};
for i=1:length(events)
    gaitEvents(i).name = eventNames{i};
    gaitEvents(i).data =events{i}';
end
%% Thoughts so ignore
% signal = LFy;invertedSignal = max(signal(:)) - signal; % get it upside donw to get mins as max
% close all
% subplot(2,1,1)
% plot(invertedSignal)
% [val loc]=findpeaks(invertedSignal); % find peak location and ind
% % hold on;plot(loc,val,'*')
% peaks(:,1) = loc;peaks(:,2) = val;
%
% %get the ones between foot strike to foot strike for each leg
% k=1;
% for i=1:length(peaks)
% if val(i) > 0.95*max(val)
%     valNew(k) = val(i);
%     locNew(k) = loc(i);
%     k=k+1;
% end
% end
% hold on;plot(locNew,valNew,'*')% get rid of those not relevant
% clear k
%
% k=1;
% SignalNew = [];% zeros(length(invertedSignal));
% SigLocNew = [];%zeros(length(invertedSignal));
% for i=1:length(invertedSignal)
% if invertedSignal(i) > 0.98*max(val)
%     SignalNew(i) = invertedSignal(i);
%     SigLocNew(i) = i;
%     k=k+1;
% else
%     SignalNew(k) = nan;
%     SigLocNew(k) = i;
%     k=k+1;
%     end
% end
% hold on;plot(SigLocNew,SignalNew);ylim([0 1000])% get rid of those not relevant
% clear k
%
% % now to get FS and TO for left foot
%
% k=1;
% for i=2:length(SignalNew)-1
%    if (isnan(SignalNew(i)) && ~isnan(SignalNew(i+1)))
%    TO(k) = SigLocNew(i);
%    k=k+1;
%    end
%
%    if (isnan(SignalNew(i)) && ~isnan(SignalNew(i-1)))
%        FS(k) = SigLocNew(i);
%   k=k+1;
%    end
%
% end
%
% hold on;plot(TO,SignalNew,'*b')
% hold on;plot(500,FS,'*g')
% % now get the gradient so that when gradient is zero get the begining and
% % end as foot strike and
% figure(2)
% c = gradient(signal)
% subplot(2,1,2)
% plot(c)
% [valGr locGr]=findpeaks(c); % find peak location and ind after gradient
% hold on;plot(locGr,valGr,'*')
%
% k=1;
% for i=1:length(valGr)
% if val(i) > 0.95*max(val)
%     valNew(k) = val(i);
%     locNew(k) = loc(i);
%     k=k+1;
% end
% end
% hold on;plot(locNew,valNew,'*')% get rid of those not relevant
% clear k
% valGr(abs(valGr) < .1)
%
% t = gradient(c);figure;plot(t)
% % d = b(a>.9*a)
% plot(b,d,'.')
% max(b)


% [LHeelStrikes,RHeelStrikes] = HeelStrikeDetection(LFy, RFy, Cutoff, SaveName)
%
% LHeelStrikesTime  = time(LHeelStrikes);
% RHeelStrikesTime = time(RHeelStrikes);
% for i=1:length(ToPlot) %
%     hold on;
%     %do this later for accuate nameing etc: [a b]=ismember('ground_force_1_vy',fields) then b =2 is fields{2}
%     h(i) = plot(ToPlot(i).data.(fields{end}),ToPlot(i).data.(fields{2}));% only the left side or do this later
%     names{i} = ToPlot(i).name;
% end
