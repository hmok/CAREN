% function gaitEvents = trcGaitEvents(trcPath)
function gaitEvents = trcGaitEvents(source)

%% Note how I have desgined the events: gaitEvents(i).data ={HSLocL, TOLocL,HSLocR, TOLocR}
% What does it do?
% it gets any trcpath (all TRCs as struct) then calcualtes gaitEvents in
% the fllowing format: Left HS, Left TO, Right HS and Right TO
% gaitEvents(i).data ={HSLocL, TOLocL,HSLocR, TOLocR};

%% TODOs
% 1. maybe I can add reproduced or DFLOW or pertubation event to this gait
% events and then I have all the possiblties of events and if I know the D
% and ND I can tell what is going on in both legs during each trial. Let;s
% think more
%% Example: run the following lines, just path is required for trc files
% addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
% source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';
% gaitEvents = trcGaitEvents(source)


% % then use it this way to plot,etc
% load gaitEvents.mat % from trcPath
% figure;title('TO Left side');vline(gaitEvents(7).data{1, 2}) % plot the TO for left side or
% figure;title('TO right side');vline(gaitEvents(7).data{1, 4}) % plot the TO for right side
% trcPath = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\trcResults';
%% main function
% x: AP in OpenSim, y=Vertical not used here, z=ML, make chagnes if the
% orierntaion of 
x=1;y=2;z=3;
trcResult = fullfile(source, 'trcResults');
cd(trcResult)
for i=1:4
    disp('be careful about the perturbed trials')
end
% cd(trcPath);
load TRCs;
for i =1:length(TRCs)
    TRCs(i).name
    sacZ = (TRCs(i).data.RPSIS(:,x)+TRCs(i).data.LPSIS(:,x))/2; %e.g. 7 means W1+L+ACC.trc
    
    % Zeni et al 2008
    % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2384115/
    
    % tHS = max(xHEEL - xSacrum)
    % tTO = min(xTOE - xSacrum)
    time = TRCs(i).data.time;
    % Left side thenright
    try
        heel2SacL = TRCs(i).data.LHEE(:,x) - sacZ;heel2SacR = TRCs(i).data.RHEE(:,x) - sacZ;
        toe2SacL = -TRCs(i).data.LMT2(:,x) + sacZ;toe2SacR = -TRCs(i).data.RMT2(:,x) + sacZ;
    catch
        if warning(['Reference to non-existent field'])
            toe2SacL = -TRCs(i).data.LMT2(:,x) + sacZ;toe2SacR = -TRCs(i).data.RMT5(:,x) + sacZ;
            
            errCollect{i} ='marker missing';
            
        end
        continue;
        
        
    end
    [~, HSLocL]= findpeaks( heel2SacL,time);[~, HSLocR]= findpeaks( heel2SacR,time);
    % figure; title(strcat('Heel 2 Sac Dis',num2str(i)));plot(heel2Sac);findpeaks(heel2Sac);
    % figure; title(strcat('Toe 2 Sac Dis',num2str(i)));plot(toe2Sac);findpeaks(toe2Sac);hold on;plot(-toe2Sac,'-')
    [~, TOLocL] = findpeaks(toe2SacL,time);[~, TOLocR] = findpeaks(toe2SacR,time);
    
    gaitEvents(i).name = TRCs(i).name;
    gaitEvents(i).data.HSLocL = HSLocL ;
    gaitEvents(i).data.TOLocL = TOLocL;
    gaitEvents(i).data.HSLocR =HSLocR;
    gaitEvents(i).data.TOLocR = TOLocR;
    gaitEvents(i).data.strideFreqL = 1./diff(HSLocL);
    %     gaitEvents(i).data.strideFreqL = 1./diff(HSLocL);
    gaitEvents(i).data.strideFreqR = 1./diff(HSLocR);
    %     gaitEvents(i).data.strideFreqR = 1./diff(HSLocR);
    
%From Motek Lua
    
%     	-- Calculate belt distance difference and distance between z of other foot at its IC and z of measured foot at its IC.
% 	local beltDifference = math.abs(beltDistance - inputs.get(sideNames[side].."Belt.Distance"))
% 	local markerDifference = math.abs(zOtherFoot - inputs.get(markerName[side]..".PosZ"))
% 	
% 	-- Step length is belt distance difference plus marker difference
% 	stepLength[side][#stepLength[side] + 1] = beltDifference + markerDifference
    
%     I do this: at IC I get diff between z of heel to heel

 gaitEvents(i).data.stepLengthL = abs((TRCs(i).data.LHEE(int64(HSLocL*100),3)-TRCs(i).data.RHEE(int64(HSLocL*100),3))/1000);
 gaitEvents(i).data.stepLengthR = abs((TRCs(i).data.RHEE(int64(HSLocR*100),3)-TRCs(i).data.LHEE(int64(HSLocR*100),3))/1000);
gaitEvents(i).data.stepWidthL = abs((TRCs(i).data.LHEE(int64(HSLocL*100),2)-TRCs(i).data.RHEE(int64(HSLocL*100),2))/1000);
 gaitEvents(i).data.stepWidthR = abs((TRCs(i).data.RHEE(int64(HSLocR*100),2)-TRCs(i).data.LHEE(int64(HSLocR*100),2))/1000);
%
%     
%     load BeltVel.mat;
%     if isa(BeltVel(i).data,'double') && length(BeltVel(i).data)>2
%         beltDisL = BeltVel(i).data(int64(HSLocL*100))*0.01;
%         beltDisR = BeltVel(i).data(int64(HSLocR*100))*0.01;
%         gaitEvents(i).data.strideLengthL = diff(TRCs(i).data.LHEE(int64(HSLocL*100),3))/1000 + beltDisL(1:end-1);
%         
%         gaitEvents(i).data.strideLengthR = diff(TRCs(i).data.RHEE(int64(HSLocR*100),3))/1000+ beltDisR(1:end-1);
%         clear beltDisL beltDisR
%     else
%         gaitEvents(i).data.strideLengthL = NaN ;
%         %     gaitEvents(i).data.strideLengthL = diff(TRCs(i).data.LHEE(HSLocL,3));
%         gaitEvents(i).data.strideLengthR = NaN;
%         %     gaitEvents(i).data.strideLengthR = diff(TRCs(i).data.LHEE(HSLocL,3));
%     end
    
    
    clear tHS HSLoc tTO TOLoc
end

save('gaitEvents','gaitEvents')


