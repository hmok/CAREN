clear
clc
%% TODOs
%1. add try catch warning if then for no pertub at all and then save them
%all as dflow.name dflow.data
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

%% main scritps
saving = 1; %save or not (1=save)
plotting = 1;% do plot
cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\trcResults')
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB2\VICON\CB2\trcResults')
load TRCs
% close all
muscelDelay=250; % ms is the muscle delay after pertubation happens
freq = 100;
slop = 3;% check this oneAcc = 0.01; %acceleation or deceleration
tdelay = 0.2; %time delay at peak
offset = 0.6; %offset from average velcoty e.g. 0.556 + offset
% Speed definition
w1=0.556; %Walk speed 2 km/s
w2=1.111; %Walk speed 4 km/s
w3=1.667; %Walk speed 6 km/s

for i=1:length(TRCs)
    
    
    try
        
        
        
        % i = 14; % which trial
        TRCs(i)
        
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
        time =  TRCs(i).data.time;
        load gaitEvents
        switch side
            case 'Left'
                a(:,1) = TRCs(i).data.LMT2(:,3)/1000; % mm 2 m
                a(:,2) = TRCs(i).data.LMT5(:,3)/1000;
                a(:,3) = TRCs(i).data.LHEE(:,3)/1000;
                HStiming = gaitEvents(i).data.HSLocL;
            case 'right'
                a(:,1) = TRCs(i).data.RMT2(:,3)/1000; % mm 2 m
                a(:,2) = TRCs(i).data.RMT5(:,3)/1000;
                a(:,3) = TRCs(i).data.RHEE(:,3)/1000;
                HStiming = gaitEvents(i).data.HSLocR;
        end
        Foot = mean(a,2);
        
        Fvel = zeros(length(Foot),1);
        Fvel = gradient(-Foot,1/freq); % left foot velocity
        % d = diff(-b)./diff(time);
        % subplot(1,2,1);
%         title(strcat('Velocity of', side, 'feet')); %plot(Fvel) %foot velcity i.e. belt velcity to reproduce the DFLOW ..
        % subplot(1,2,2); title('Acceleration of feet');plot(c)
        
        %% Let's reproduce Belt velocity from foot markers and velcity
        % i have the peak so that is where velcoty of feet must match and then
        % recreate the following pattern then I have the timing of pertub
        % initiation, then I can reproduce the next graph (as i have the length of time e.g. length(time) =700 or any)
        
        % find peak of velcity
        % [LmaxVel Lloc]=max(-Fvel);
%         p=1;
%         for l=1:2:length(HStiming)-1
%         [Pks,LocS] = findpeaks(abs(Fvel(HStiming(l)*100:HStiming(l+1)*100)))    
%         LocSS(:,p) = LocS
%         p=p+1;%
%         end
        [Pks,LocS] = findpeaks(abs(Fvel));
        [~, Location] = min(abs(abs(Fvel(LocS))-(vel+offset)));
        Lloc = LocS(Location);
        
        if plotting == 1 && ~strcmp(side, 'No') %#ok<STCMP>
            figure(i);
            findpeaks(abs(Fvel));
        end
        
        tPeak = time(Lloc); % time of the peak vecolity
        delta = offset/slop;  % Acc = tan(theta) = 0.01 = offset/delta(t)
        % we need 6 points (X =time,Y=velocity) to reproduce dflow belt speed
        frame3 = Lloc - 0.5*tdelay*freq; %- muscelDelay; % this is the frame # of third point
        frame4 = Lloc + 0.5*tdelay*freq; % - muscelDelay; % this is the frame # of third point
        frame2 = frame3-delta*freq;
        frame5 = frame4+delta*freq;
        frame6 = length(Fvel);
        % vline(frame4) (frame4-frame3)/freq
        dflow(1:frame2) = vel;
        
        %points between frame2 and 3
        NumberNewPoints = frame3-frame2-2;
        xvals = round(linspace(frame2, frame3, NumberNewPoints+2));
        yvals = linspace(vel,vel+offset, NumberNewPoints+2);
        pts = [xvals(:), yvals(:)];
        dflow(frame2:frame3-1) = pts(:,2);
        
        
        %points between frame4 and 5
        NumberNewPoints = frame5-frame4-2;
        xvals = round(linspace(frame4, frame5, NumberNewPoints+2));
        yvals = linspace(vel+offset,vel, NumberNewPoints+2);
        pts = [xvals(:), yvals(:)];
        dflow(frame4:frame5-1) = pts(:,2);
        
        dflow(frame3:frame4) = vel+offset;%check or LmaxVel = w1+offset or LmaxVel
        dflow(dflow<0) = 0;%check or LmaxVel = w1+offset or LmaxVel
        dflow(frame5:frame6) = vel; %change LmaxVel = w1+offset
        
        
        
        % ws = [w1 w1 w1+offset w1+offset w1 w1];
        % frames = [1 frame2 frame3 frame4 frame5 frame6];
        %
        %
        % for i=1:6
        % hold on;plot(frames(i),ws(i),'.-')
        % end
        % plot(frames,ws)
        if plotting == 1 && ~strcmp(side, 'No')
            hold on; plot(dflow)
        end
    catch
        if warning(['This is normal walking or not Pertub, check the name and if doubt check the c3d file '])
%             errCollect(i).name = TRCs(i).name;% c3dp(k).name; %Store the files names with error it also finalizes the script
%             errCollect(i).data = NaN;% c3dp(k).name; %Store the files names with error it also finalizes the script
%             k=k+1;
            BeltVel(:,i).name = TRCs(i).name;
            BeltVel(:,i).data = NaN;
            
        end
        continue;
    end
    BeltVel(:,i).name = TRCs(i).name;
    BeltVel(:,i).data = dflow';
    
end

if saving == 1
save('BeltVel', 'BeltVel')
end