%% Lets do Margine of Stability (MoS)

% function [BoS, MoS, SL, SF, CoM] = gaitUtils(trcPath)
%% TODOS
% 1.there are several issues: here get the equation corrected. 
% equations check wi lit
% graph plotted and data saved make sure BOS, SL, SF, CoM, etc are being
% calculated...
% save it as follows all the output in gaitPar = {mos_AP, mos_ML, bos_AP,
% bos_ML, com, stLength, stFreq}; how?
% 1. input: trc files, locations

% 2. outputs: get COM, SL, MOS, BS, ...
% 3. 

% Author(s): Hossein Mokhtarzadeh
% Contibutor(s): Carlos Eduardo Landeau Bobadilla 

%% path etc
trcPath = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\trcResults';
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
cd(trcPath);
load TRCs;
%% get freq of one c3d in the list of c3ds

[ParentFolderPath] = fileparts(trcPath);
[~, ParentFolderName] = fileparts(ParentFolderPath);
c3dpath = fullfile(ParentFolderPath,'C3DFiles','W101.c3d'); % this may need to change if naming changes

% Construct an opensimC3D object with input c3d path
% Constructor takes full path to c3d file and an integer for forceplate
% representation (1 = COP). 

c3d = osimC3D(c3dpath,1);
% Get some stats...
% Get the marker data rate
rMarkers = c3d.getRate_marker();
freq = rMarkers;

%% direction in OpenSim based of CAREN
% clear 
x=1;% x: Left Medial (M) (1)
y=2;% y: Vertical Up (2)
z=3;% z: Forward Anterior (3)
% freq=100;


%% loop to get CoM


for i =11%:9%1:length(TRCs)
    
    disp('this is being processed:')
    TRCs(i).name
%     sacZ = (TRCs(i).data.RPSIS(:,3)+TRCs(i).data.LPSIS(:,3))/4; %e.g. 7 means W1+L+ACC.trc
com1 = (TRCs(i).data.LASIS/1000+TRCs(i).data.RASIS/1000+TRCs(i).data.LPSIS/1000+TRCs(i).data.RPSIS/1000); %mm to m
com1 = com1/4;
% check whether com is 700 or similar to Belt Velocty ones
% C_O_M = (C_o_M_matrix)/4; %Centre of mass 

    
    % Zeni et al 2008
    % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2384115/
    
    % tHS = max(xHEEL - xSacrum)
    % tTO = min(xTOE - xSacrum)
    time = round(TRCs(i).data.time(end));
    t = TRCs(i).data.time;
% I am correcting any missign data in TRCs as some of them are 699 or 698
% frames, good idea to correct it when c3d are being converted. check that
% later
xx = 0:1/freq:(time-1/freq);
for p=1:3
comX(:,p) = spline(t,com1(:,p),xx);
end
com = comX; % now the COM is in the same no as other file like DFLOW.


    
    % Left side thenright
    
    %% Calculte MoS 



%% Calculate pendulum's length and determine eigenfrequency

l = max(com(:,y));% y is vertical or 2
inv_pend_eigenfreq = sqrt(l / 9.81); % g =9.81 m/s^2 

%% Calculate XCoM (for x and z directions) need thought check 
% Get the velocity from TRC, check the equations all the details

com_vel(:,x) = gradient(-com(:,x),1/freq);
com_vel(:,y) = gradient(-com(:,y),1/freq);
com_vel(:,z) = gradient(-com(:,z),1/freq);
plot(com_vel(:,x))

% com_vel_x = diff(com(:, 1)/1000)./diff(markerStruct.(fields{end}));
xcom_x = com(:, x) + (com_vel(:,x) * inv_pend_eigenfreq);
xcom_y = com(:, y) + (com_vel(:,y) * inv_pend_eigenfreq);
xcom_z = com(:, z) + (com_vel(:,z) * inv_pend_eigenfreq);

% com_z_tm = com(:, 3) * -1; 
% com_vel_z_tm = diff(com(:, 3)/1000)./diff(markerStruct.(fields{end}));%differentiate_data(com(:, 3));
%  xcom_z_tm = com_z_tm(1:end-1, 1)/1000 + (com_vel_z_tm * inv_pend_eigenfreq);

% Belt speed from TRC or DFLOw, let;s do TRC only for now and then I will
% compare them with DFLOW.
load BeltVel.mat

if contains(BeltVel(i).name,TRCs(i).name,'IgnoreCase',true)
    belt_velz = BeltVel(i).data;
end
belt_dist =  belt_velz;
if size(com,1) == size(belt_dist,1) 
com(:, 3) = com(:, 3) * -1 + belt_dist; % convert to "overground" values
elseif size(com,1) > size(belt_dist,1) 
    com(1:end-1, 3) = com(1:end-1, 3) * -1 + belt_dist; % convert to "overground" values
end

%% z: forward (AP), x: to left (ML), y: up (vertical) as in OpenSim careful about z and x which are diff from openSim vs. CAREN
% Calcualte MOS, BOS (based of Supportt)
mos_AP = []; %zeros(1,length(TRCs(i).data.time));
mos_ML = []; %zeros(1,length(TRCs(i).data.time));
% for i = 1:length(TRCs(i).data.LASIS)
bos_AP_L = TRCs(i).data.LMT2(:,z)/1000 ;bos_ML_L = TRCs(i).data.LLM(:,x)/1000;
bos_AP_R = TRCs(i).data.RMT2(:,z)/1000 ;bos_ML_R = TRCs(i).data.RLM(:,x)/1000;

if bos_AP_L > bos_AP_R % TRCs(i).data.LMT2(i,3) > TRCs(i).data.RMT2(i,3) % MoS
    mos_AP(:,i) = bos_AP_L - xcom_z; %TRCs(i).data.LMT2(:,3)/1000 - xcom_z; 
    mos_ML(:,i) = bos_ML_L - xcom_x; % TRCs(i).data.LLM(i,3)/1000 - xcom_z;
    bos_AP(:,i) = bos_AP_L;
    bos_ML(:,i) = bos_ML_L;
else
    mos_AP(:,i) = bos_AP_R - xcom_z; %TRCs(i).data.RMT2(:,3)/1000 - xcom_z;   
    mos_ML(:,i) = bos_ML_R - xcom_x; %TRCs(i).data.RLM(:,3)/1000 - xcom_z;   
     bos_AP(:,i) = bos_AP_R;
    bos_ML(:,i) = bos_ML_R;
end 
% end
    
%     gaitEvents(i).name = TRCs(i).name;
%     gaitEvents(i).data ={HSLocL, TOLocL,HSLocR, TOLocR};
%     clear tHS HSLoc tTO TOLoc
% gaitPar = {mos_AP, mos_ML, bos_AP, bos_ML, com, strideLength, strideFreq};

% to get stride length and freq
load gaitEvents % get all the trc files etc
% strideLengthL =   % distance each heel travles and its 1/time i.e. freq
% strideFreqL =  % distance each heel travles and its 1/time i.e. freq
% strideLengthR =  % distance each heel travles and its 1/time i.e. freq
% strideLengthL =  % distance each heel travles and its 1/time i.e. freq

    gaitPar(i).name = TRCs(i).name;
    gaitPar(i).data.mos_AP = mos_AP ; 
    gaitPar(i).data.mos_ML = mos_ML;
    gaitPar(i).data.bos_AP =bos_AP;
    gaitPar(i).data.bos_ML = bos_ML;
    gaitPar(i).data.com = com;
    gaitPar(i).data.strideLength = strideLength;
    gaitPar(i).data.strideFreq = strideFreq;
    gaitPar(i).data.bos_ML = bos_ML;
    clear mos_AP mos_ML bos_AP bos_ML com strideLength strideFreq
end

%% saving all in gaitPar
% gaitPar = {mos_AP, mos_ML, bos_AP, bos_ML, com, stLength, stFreq};

save('gaitPar','gaitPar')

%% plotting

%% plotting mos
limit = [-1 1];xlimi = [0 time];
subplot(2,2,1); plot(TRCs(i).data.time,mos_AP);title('AP MoS');ylim(limit);xlabel('Time (s)');xlim(xlimi)
subplot(2,2,2); plot(TRCs(i).data.time,mos_ML);title('ML MoS');ylim(limit);xlabel('Time (s)');xlim(xlimi)
subplot(2,2,3); plot(TRCs(i).data.time,bos_AP);title('AP BoS');ylim(limit);xlabel('Time (s)');xlim(xlimi)
subplot(2,2,4); plot(TRCs(i).data.time,bos_ML);title('ML BoS');ylim(limit);xlabel('Time (s)');xlim(xlimi)

load gaitEvents

figure;
hold on;plot(TRCs(i).data.time,bos_AP);title('AP MoS'); %ylim(limit);xlabel('Time (s)');xlim(xlimi)
vline(gaitEvents(i).data.HSLocR(:)); %title('AP MoS'); %ylim(limit);xlabel('Time (s)');xlim(xlimi)

[~,tt]=max(mos_AP(int(gaitEvents(i).data.HSLocR*100),10))
mos_AP(int64(gaitEvents(i).data.HSLocR*100),11)
% subplot(2,1,1);hold on;plot(mos_AP);legend('mos AP')
% subplot(2,1,2);hold on;plot(mos_ML);legend('mos ML')


% save('gaitEvents','gaitEvents')


