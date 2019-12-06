
%  function dflowReproduceTRC_Events()
%% TODOs
%1. it is quite ok but the algorithm is not perfect yet
% ask the user for eaech if corect let go if not correct it and then do the
% rest for frame and save it...ayval
% add try catch warning if then for no pertub at all and then save them
%all as dflow.name dflow.data
%add to correct them and save them all with the names and that is it...
%nans for no pertub, compr with actual data
% error to be added compare with name of the file...
% add 26 27 CB1 wrong I think the slop offset must be correctd
% this is only for Left ACC so make it for all possiblities L R and ACC DEC
%try and catch for w1 +ACC or DEC
%check errors which are mostly realted to DEC esp. W3 and W2

% States & Variables:

%              |  Time  |
%     -        ..........
%             .          .
%   offset   .            .
%           .              .
% ..........                ...........
% | 0  | 1 | 2 |   3    | 4 |   0       These are the different states per DoF / Belt


% Slope is determined by Acc

function [predictedFrames, BeltVelPredict,FootVelL, FootVelR] = dflowReproduceTRC_alignSignal(trcPath, saving, plotting)


%% note that BeltVel can be slightly different from what actually happedn but it is a good approximation. This can be later
% imporved. so one can use this to input trcPath and get timing of when
% pertubation occured, the frame number and overal BeltVel in all trc files
% in the folder input.
%% Example:
% [framePerturb BeltVel] = dflowReproduceTRC_Events(trcPath, saving (1 or -1), plotting (1, -1))
% trcPath = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\DF1\trcResults';
% trcPath = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\trcResults';
% 
% 1. addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab')
% 2. trcPath = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\CB1\trcResults';

% CHOOSE between 3,4,5 below

%  3. [framePerturb, BeltVelPredict,FootVelL, FootVelR] = dflowReproduceTRC_alignSignal(trcPath, -1, -1);
%or 4. [framePerturb, BeltVelPredict,FootVelL, FootVelR] = dflowReproduceTRC_alignSignal(trcPath, -1, 1); %plot Yes


%or CAREFul SAVING AS WELL: 
% 5. [framePerturb, BeltVelPredict,FootVelL, FootVelR] = dflowReproduceTRC_alignSignal(trcPath,1, 1); %plot, SAVE
% 

%% main scritps
% saving = 1; %save or not (1=save)
% plotting = 1;% do plot
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\trcResults')
cd(trcPath)
load TRCs
close all
cd ../ 
%dflowDir = strcat(cd,'\DFLOW'); 
dflowDir = cd;
disp(dflowDir);

cd('C:\Users\Corrine\Desktop\Intern2019Summer\OpenSim\Matlab\mFiles')
[~,~, PerturbFrame] = dflow(dflowDir, -1, -1); % this is from actual data
%[~,~, PerturbFrame] = dflow(dflowDir, -1); % this is from actual data
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB2\VICON\CB2\trcResults')
% load gaitEvents.mat
% close all
muscelDelay=.30; % 200-250ms is the muscle delay after pertubation happens so it becomes .2(s)*freq(100)=20
freq = 100;
slop = 3;% check this oneAcc = 0.01; %acceleation or deceleration
tdelay = 0.2; %time delay at peak
offset = 0.6; %offset from average velcoty e.g. 0.556 + offset
% Speed definition
w1=0.556; %Walk speed 2 km/s
w2=1.111; %Walk speed 4 km/s
w3=1.667; %Walk speed 6 km/s
vel = 0;
for i=1:length(TRCs)
    
    
    %     try
    
    
    
    % i = 14; % which trial
    TRCs(i)
    i
    framePerturb(:,i).name = TRCs(i).name;
    %% W1 and L R ACC DEC and W2, W3
    Ws = {'W1','W2','W3'};
    Vs = [w1, w2, w3];
    k=1;
    side = 'No';
    for k=1:3 % velocties
        if contains(TRCs(i).name,Ws(k)) && contains(TRCs(i).name,'+ACC') && contains(TRCs(i).name,'+R')
            vel = Vs(k); slop = 3; side = 'right'; offset = 0.6;
            disp(strcat('Vel: ', num2str(vel),' Slope: ',num2str(slop),' Side: ',side))
            %                 Perturb = 1; % there is a pertub so plot it if asked
        elseif contains(TRCs(i).name,Ws(k)) && contains(TRCs(i).name,'+DEC') && contains(TRCs(i).name,'+R')
            vel = Vs(k); slop = -3; side = 'right';offset = -0.6;
            disp(strcat('Vel: ', num2str(vel),' Slope: ',num2str(slop),' Side: ',side))
            %                 Perturb = 1;
        elseif contains(TRCs(i).name,Ws(k)) && contains(TRCs(i).name,'+ACC') && contains(TRCs(i).name,'+L')
            vel = Vs(k); slop = 3; side = 'Left';offset = 0.6;
            disp(strcat('Vel: ', num2str(vel),' Slope: ',num2str(slop),' Side: ',side))
            %                 Perturb = 1;
        elseif contains(TRCs(i).name,Ws(k)) && contains(TRCs(i).name,'+DEC') && contains(TRCs(i).name,'+L')
            vel = Vs(k); slop = -3; side = 'Left';offset = -0.6;
            disp(strcat('Vel: ', num2str(vel),' Slope: ',num2str(slop),' Side: ',side))
            %                 Perturb = 1;
        else
            disp(strcat('Not in Vel: ', Ws{k}));
            disp(TRCs(i).name);
            
        end
    end
    
    %%
    disp(side)
    time =  TRCs(i).data.time;fnames=fields(TRCs(i).data);
    HStiming =[];
    %         load gaitEvents.mat; % to find intersection in between left right
    %         switch side
    %             case 'Left'
    a(:,1) = TRCs(i).data.LMT2(:,3)/1000; % mm 2 m
    a(:,2) = TRCs(i).data.LMT5(:,3)/1000;
    a(:,3) = TRCs(i).data.LHEE(:,3)/1000;
    %                 HStiming=gaitEvents(i).data.HSLocL; % this is the HS timing for left foot
    %                 stepL=gaitEvents(i).data.stepLengthL;
    %     try
    if contains('RMT2', fnames)
        b(:,1) = TRCs(i).data.RMT2(:,3)/1000; % mm 2 m
    else
        b(:,1) = TRCs(i).data.RMT5(:,3)/1000; % mm 2 m
    end
    b(:,2) = TRCs(i).data.RMT5(:,3)/1000;
    b(:,3) = TRCs(i).data.RHEE(:,3)/1000;
    %                  Events=gaitEvents(i).data.HSLocL*freq; % this will provide where pertub or gaite events hapn to find >2 intersection that is pertub
    %               case 'right'
    a(:,1) = TRCs(i).data.LMT2(:,3)/1000; % mm 2 m
    a(:,2) = TRCs(i).data.LMT5(:,3)/1000;
    a(:,3) = TRCs(i).data.LHEE(:,3)/1000;
    
    
    if contains('RMT2', fnames)
        b(:,1) = TRCs(i).data.RMT2(:,3)/1000; % mm 2 m
    else
        b(:,1) = TRCs(i).data.RMT5(:,3)/1000; % mm 2 m
    end
    %     b(:,1) = TRCs(i).data.RMT2(:,3)/1000; % mm 2 m
    b(:,2) = TRCs(i).data.RMT5(:,3)/1000;
    b(:,3) = TRCs(i).data.RHEE(:,3)/1000;
    %                 HStiming=gaitEvents(i).data.HSLocR; % this is the HS timing for right foot
    %                 stepL=gaitEvents(i).data.stepLengthR;
    %                  Events=gaitEvents(i).data.HSLocR*freq;
    %
    %         end
    FootL = mean(a,2);
    FootR = mean(b,2);
    clear a b
    FvelL = gradient(-FootL,1/freq); % left foot velocity
    FvelR = gradient(-FootR,1/freq); % left foot velocity
    if plotting ==1
    figure; plot(FvelL);hold on;plot(FvelR);hold off;legend('Left','Right')
    end
    % hold on;vline(gaitEvents(i).data.HSLocL*freq)
    % hold on;vline(gaitEvents(i).data.HSLocR*freq)
    % hold on;plot(FvelR);plot(smooth(FvelR, 2, 100, 4)); legend('raw','filtered')
    % ipt = findchangepts(FvelL)
    
    %     FvelR=smooth(FvelR, 4, 100, 4);
    %     FvelL=smooth(FvelL, 4, 100, 4);
    %      figure; plot(FvelL);hold on;plot(FvelR);hold off;legend('Left','Right')
    X1=[1:length(FvelL)]';X2=X1;pertubL=[];pertubR=[];frame2=[];
    
    %     hold on; plot(sqrt(FvelL.^2)); plot(sqrt(FvelR.^2));legend('Left','Right');hold off
    cd('C:\Users\Corrine\Desktop\Intern2019Summer\OpenSim\Matlab\mFiles');
    intSec = intersections(X1,FvelL,X2,FvelR);
    
    % a1=FvelL;
    %  plot(a1)
    %  a1_fft = abs(fft(a1));
    % a1_fft = a1_fft(1:length(a1)/2);
    % plot(a1_fft);
    % X = ifft(a1_fft);plot(X)
    % vline(intSec)
    %     find abs of max of both signals L R
    %         maxFL = max(abs(FvelL));maxFR = max(abs(FvelR));
    
%     thresholdDEC=.95;threshold2=1.6;thresholdACC=.5;
%     for j=1:length(intSec)-1
%         % numberof sample btwn events >30 and
%         h= round(intSec(j):intSec(j+1));
%         %         if length(h)< 2*mean(diff(intSec))
%         
%         %         if strcmp(side,'Left') && ((abs(mean(FvelL(h))) - vel)>offset/2) && abs(abs(mean(FvelR(h))) - vel)>.2*vel
%         if strcmp(side,'Left') && contains(TRCs(i).name,'+ACC') && max((abs(FvelL(h))) - vel)>thresholdACC*abs(offset)  && max(abs(FvelR(h)))>threshold2*vel
%             
%             frame2 = round(intSec(j));
%             break
%             %         elseif strcmp(side,'right') && ((abs(mean(FvelR(h))) - vel)>offset/2) &&  abs(abs(mean(FvelL(h))) - vel)>.2*vel
%             
%         elseif strcmp(side,'Left') && contains(TRCs(i).name,'+DEC') && ( abs(min((abs(FvelL(h))) - vel))>thresholdDEC*abs(offset)) && max(abs(FvelR(h)))>threshold2*vel
%             frame2 = round(intSec(j));
%             break
%             
%         elseif strcmp(side,'right') && contains(TRCs(i).name,'+ACC') && (max(abs(FvelR(h)) - vel)>thresholdACC*abs(offset) ) && max(abs(FvelL(h)))>threshold2*vel
%             frame2 = round(intSec(j));
%             break
%             
%         elseif strcmp(side,'right') && contains(TRCs(i).name,'+DEC') && abs(min((abs(FvelR(h))) - vel))>thresholdDEC*abs(offset)  && max(abs(FvelL(h)))>threshold2*vel
%             frame2 = round(intSec(j));
%             break
%         end
%         
%         clear h
%         %         end
%         
%     end
    
    %%     New idea
    % % 1. to get everthing minus vel to start from zero
    % 2. make abs of all
    % 3. get all below (offset + .1) and and its max, this is the peak
    % 4. Pertub is the first closest before the time of max in #3
    
    FvelLabs = abs(FvelL - vel);FvelRabs = abs(FvelR - vel); %1,2
    % figure(2);hold on; plot(FvelLabs); plot(FvelRabs);legend('ABS Left','ABS Right');
    %  o=FvelRabs(FvelRabs>.5);findpeaks(o)
    if strcmp(side,'Left')
        %     [PKS,LOCS]=findpeaks(smooth(FvelLabs,3,100,4));
        [PKS,LOCS]=findpeaks(FvelLabs);
        [minpk,~]= max(PKS(PKS<(abs(offset)+.1) & PKS>(abs(offset)-.2) ));%2*abs(offset))); %smooth(FvelR, 4, 100, 4)
        [lia,locb]=ismember(minpk,PKS);
        if ~isempty(locb)
            jj=find(intSec>LOCS(locb),1)-1;
            frame2 =round(intSec(jj));
        end
        %     for h=1:100
        %         if abs(FvelL(LOCS(locb)-h)-vel) < .1
        %             frame2 = LOCS(locb)-h;
        %             break
        %         end
        %
        %     end
        clear PKS LOCS lia locb h minpk jj
    elseif strcmp(side,'right')
        %     [PKS,LOCS]=findpeaks(smooth(FvelRabs,3,100,4));
        [PKS,LOCS]=findpeaks(FvelRabs);
        [minpk,~]= max(PKS(PKS<(abs(offset)+.1) & PKS>(abs(offset)-.2) ));%2*abs(offset))); %smooth(FvelR, 4, 100, 4)
        [lia,locb]=ismember(minpk,PKS);
        %     for h=1:100
        %         if abs(FvelR(LOCS(locb)-h)-vel) < .1
        %             frame2 = LOCS(locb)-h;
        %             break
        %         end
        %
        %     end
        if ~isempty(locb)
            jj=find(intSec>LOCS(locb),1)-1;
            
            frame2 =round(intSec(jj));
        end
        clear PKS LOCS lia locb h minpk jj
    end
    
    %
    tt = strcat(side,' : ', int2str(frame2),TRCs(i).name);title(tt);
    if ~isempty(frame2) && plotting ==1
        vline(frame2)
    else
        sttr = strcat('there was no pertubation in this trial: ',TRCs(i).name);
        disp(sttr)
        erroCollect(:,i).name = extractBefore(TRCs(i).name,'.');
        
    end
    framePerturb(:,i).name = extractBefore(TRCs(i).name,'.');
    
    %% Perturbation prediction and actual one
    if ~strcmp(side, 'No')
    
          % bring the pertub fro actual data here as well and diplay it
          disp(length(PerturbFrame));
          for jj = 1:length(PerturbFrame)
              
              index = find(strcmp({framePerturb(:,i).name}, {PerturbFrame(jj).name})==1); % PerturbFrame is the actual one frmom dflow.m
              if ~isempty(index)
                  disp(strcat('This is the actual frame at Pertubation: ', int2str(PerturbFrame(jj).iniPerturb)));
                  ActualFrame2 = PerturbFrame(jj).iniPerturb;
                  break
              else
                  ActualFrame2 = [];
              
              end
              
          end
          prompt = 'Is this correct i.e. Pertubation has been predicted correctly (Y=1/NA=3)? ';
          x = input(prompt);
        
      %%  these are options for manual and automatic handingling of actual vs. pridicted pertubation timing
%       x = 1 means either predicted or actual exist so the actual one will be presented if the actual one is not there so we ask for prediction so input does it
%         x=3 means the pertub did not happn or they don;t match at allthe
%         i.e.       actual vs. predicted.
if (x == 1) && isempty(ActualFrame2)
    disp('find it manually and then')
    frame2 = input('Type the actual frame location i.e. 213:  ');
    framePerturb(:,i).data = frame2; %actual from
    disp(~isempty(ActualFrame2))
    if x == 1 && ~isempty(ActualFrame2)
        framePerturb(:,i).data = ActualFrame2; %actual from
        frame2 = ActualFrame2;
    end
elseif x == 3
    framePerturb(:,i).data = NaN;
    
else
    framePerturb(:,i).data = frame2;
end
    
     end
    
    
    %         if plotting == 1 && ~strcmp(side, 'No') %#ok<STCMP>
    %             figure(i);
    %             findpeaks(abs(Fvel));
    %         end
    
    
    %why 7 ???????
   
    dflowPredic = vel*ones(length(TRCs(1).data.time),1);
    if ~isempty(frame2)
        %         tPeak = time(Lloc); % time of the peak vecolity
        delta = offset/slop;  % Acc = tan(theta) = 0.01 = offset/delta(t)
        % we need 6 points (X =time,Y=velocity) to reproduce dflow belt speed
        frame3 = frame2 + delta*freq; %tdelay*freq - muscelDelay; % this is the frame # of third point
        frame4 = frame3 + tdelay*freq; % this is the frame # of third point
        %         frame2 = frame3-delta*freq;
        frame5 = frame4+delta*freq;
        frame6 = length(FvelR);
        % vline(frame4) (frame4-frame3)/freq
        dflowPredic(1:frame2) = vel;
        
        %points between frame2 and 3
        NumberNewPoints = frame3-frame2-2;
        xvals = round(linspace(frame2, frame3, NumberNewPoints+2));
        yvals = linspace(vel,vel+offset, NumberNewPoints+2);
        pts = [xvals(:), yvals(:)];
        dflowPredic(frame2:frame3-1) = pts(:,2);
        
        
        %points between frame4 and 5
        NumberNewPoints = frame5-frame4-2;
        xvals = round(linspace(frame4, frame5, NumberNewPoints+2));
        yvals = linspace(vel+offset,vel, NumberNewPoints+2);
        pts = [xvals(:), yvals(:)];
        dflowPredic(frame4:frame5-1) = pts(:,2);
        
        dflowPredic(frame3:frame4) = vel+offset;%check or LmaxVel = w1+offset or LmaxVel
        dflowPredic(dflowPredic<0) = 0;%check or LmaxVel = w1+offset or LmaxVel
        dflowPredic(frame5:frame6) = vel; %change LmaxVel = w1+offset
        
    end
    
    
    % ws = [w1 w1 w1+offset w1+offset w1 w1];
    % frames = [1 frame2 frame3 frame4 frame5 frame6];
    %
    %
    % for i=1:6
    % hold on;plot(frames(i),ws(i),'.-')
    % end
    % plot(frames,ws)
    if plotting == 1 && ~strcmp(side, 'No')
        hold on; plot(dflowPredic)
    end
    %     catch
    %         if warning(['This is normal walking or not Pertub, check the name and if doubt check the c3d file '])
    % %             errCollect(i).name = TRCs(i).name;% c3dp(k).name; %Store the files names with error it also finalizes the script
    % %             errCollect(i).data = NaN;% c3dp(k).name; %Store the files names with error it also finalizes the script
    % %             k=k+1;
    % %             BeltVel(:,i).name = TRCs(i).name;
    % %             BeltVel(:,i).data = NaN;
    %
    %         end
    %         continue;
    %     end
    BeltVel(:,i).name = extractBefore(TRCs(i).name,'.');
    BeltVel(:,i).data = dflowPredic';
    FootVelR(:,i).name = extractBefore(TRCs(i).name,'.');
    FootVelR(:,i).data = FvelR;
    FootVelL(:,i).name = extractBefore(TRCs(i).name,'.');
    FootVelL(:,i).data = FvelL;
    clear dflowPredic
end
% saving

BeltVelPredict = BeltVel;
predictedFrames = framePerturb;
if saving == 1
cd(trcPath)
    save('BeltVelPredict', 'BeltVelPredict')
    save('predictedFrames', 'predictedFrames')
end