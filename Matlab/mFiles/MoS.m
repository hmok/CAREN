%% Lets do Margine of Stability (MoS)

% function [BoS, MoS,CoM, XCoM] = gaitUtils(trcPath)
function [mos_AP mos_ML bos_AP bos_ML, com_vel XCOM COM] = MoS(source,saving)

%% TODOS
% there are several issues: here get the equation corrected.
% 1. check the equtions and based them on sth solid a paper..and compare
% then check the Gordon 2018 and Hof et al 2005, The condition for dynamic
% stability read and get the details perfectly and code them here nicely
% 2. Check all the details I do correct for elemtn size but still have probelm alter so equation and etc
% input: trc files, locations
% 2. outputs: get COM, SL, MOS, BS, ...
% 3.

% Author(s): Hossein Mokhtarzadeh
% Contibutor(s): Carlos Eduardo Landeau Bobadilla

%% Example:
% 1. provide the source folder
% source ='C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';
% 2. run this beaut below:
%  [mos_AP mos_ML bos_AP bos_ML, com_vel XCOM COM] = MoS(source,-1)
%% direction in OpenSim based of CAREN
%  clear
x=1;% x: Left Medial (M) (1)
y=2;% y: Vertical Up (2)
z=3;% z: Forward Anterior (3)
freq=100;

%% path etc and belt velocity
% source ='C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';
% trcPath = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\trcResults';
trcPath = fullfile(source,'trcResults');
cd(trcPath);
load TRCs;
% load('C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\DFLOW\dflow.mat')
if exist(fullfile(source,'DFLOW','dflow.mat'))
    load(fullfile(source,'DFLOW','dflow.mat'));
else
    disp('This source has no DFLOW data, correct it and get back again: ')
    source
    
end

%% loop to get CoM

% for i =6:length(TRCs)     why 6 ???????
for i =2:length(TRCs)% ignore CALIBRATION!!!!!! different dimension£¡
    
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
    
    
    [~,name,~]=fileparts(TRCs(i).name);
    
    %% Belt Speed
    disp(length(dflow.BeltVel))
    for j=1:length(dflow.BeltVel)
        
        index1 = find(strcmp({dflow.BeltVel.name}, name)==1);
        % index2 = find(strcmp({BeltVelPredict.name}, trials(i))==1);
        
        if  ~isempty(index1)
            belt_velz = dflow.BeltVel(index1).data;
        else
            disp('This is missing check this out')
            %     W = {'W1','W2','W3'};
            disp(name)
            switch name
                
                case contains(name,'W1')
                    belt_velz = 2/3.6*ones(length(comX),1);
                    disp(belt_velz)
                case contains(name,'W2')
                    belt_velz = 4/3.6*ones(length(comX),1);
                    disp(belt_velz)
%                 dsp vs disp???????????
                case contains(name,'W3')
                    belt_velz = 6/3.6*ones(length(comX),1);
                    disp(belt_velz)
                otherwise %deal with calibration
                     belt_velz = 0;
                     disp(belt_velz)
            end
        end
    end
    
    %% Calculate XCoM (for x and z directions)
    % Get the velocity from TRC, check the equations all the details
    
    com_vel(:,x) = gradient(-com(:,x),1/freq);
    com_vel(:,y) = gradient(-com(:,y),1/freq);
    com_vel(:,z) = gradient(-com(:,z),1/freq);
    % plot(com_vel(:,x))
    COM_VEL(:,:,i) = com_vel;
    % com_vel_x = diff(com(:, 1)/1000)./diff(markerStruct.(fields{end}));
    
    
    xcom_x = com(:, x) + (com_vel(:,x) * inv_pend_eigenfreq);
    xcom_y = com(:, y) + (com_vel(:,y) * inv_pend_eigenfreq);
    xcom_z = com(:, z) + (com_vel(:,z) * inv_pend_eigenfreq);
    
    XCOM(:,1,i) = xcom_x;
    XCOM(:,2,i) = xcom_y;
    XCOM(:,3,i) = xcom_z;
    
    COM(:,1,i) = com(:, x);
    COM(:,2,i) = com(:, y);
    COM(:,3,i) = com(:, z);
    
    % COM_x(:,i) = xcom_x;
    % com_z_tm = com(:, 3) * -1;
    % com_vel_z_tm = diff(com(:, 3)/1000)./diff(markerStruct.(fields{end}));%differentiate_data(com(:, 3));
    %  xcom_z_tm = com_z_tm(1:end-1, 1)/1000 + (com_vel_z_tm * inv_pend_eigenfreq);
    
    % Belt speed from TRC or DFLOw, let;s do TRC only for now and then I will
    % compare them with DFLOW.
    % load BeltVel.mat
    
    
    % trials = {'W1+L+ACC';'W1+L+DEC';'W1+R+ACC';'W1+R+DEC';'W1+T+L+ACC';'W1+T+L+DEC';'W1+T+R+ACC';'W1+T+R+DEC';...
    %     'W2+L+ACC';'W2+L+DEC';'W2+R+ACC';'W2+R+DEC';'W2+T+L+ACC';'W2+T+L+DEC';'W2+T+R+ACC';'W2+T+R+DEC';...
    %    'W3+L+ACC';'W3+L+DEC';'W3+R+ACC';'W3+R+DEC';'W3+T+L+ACC';'W3+T+L+DEC';'W3+T+R+ACC';'W3+T+R+DEC' };
    % belt_dist =  belt_velz;
    % if size(com,1) == size(belt_velz,1)
    
    % elseif size(com,1) > size(belt_velz)
    %     com(1:end-1, 3) = com(1:end-1, 3) * -1 + belt_velz; % convert to "overground" values
    % end
    
    com(:, 3) = com(:, 3) * -1 +  cumtrapz(belt_velz)*1/freq; %belt_velz; % convert to "overground" values beltDis = cumtrapz(belt_velz(1:699,:))*1/freq;
    
    %% z: forward, x: to left, y: up
    % Calcualte MOS
    % mos_AP = zeros(1,length(markerStruct.(fields{end})(1:end-1)));
    % mos_ML = zeros(1,length(markerStruct.(fields{end})(1:end-1)));
    % for i = 1:length(TRCs(i).data.LASIS)
    
    LMT2AP = TRCs(i).data.LMT2(:,3)/1000; %AP
    LLMAP = TRCs(i).data.LLM(:,3)/1000;
    RMT2AP = TRCs(i).data.RMT2(:,3)/1000;
    RLMAP = TRCs(i).data.RLM(:,3)/1000;
    LMT2ML = TRCs(i).data.LMT2(:,1)/1000; %AP
    LLMML = TRCs(i).data.LLM(:,1)/1000;
    RMT2ML = TRCs(i).data.RMT2(:,1)/1000;
    RLMML = TRCs(i).data.RLM(:,1)/1000;
    if length(xcom_z) ~= length(TRCs(i).data.LMT2(:,3))
        LMT2AP = interp1(1:numel(LMT2AP), LMT2AP, linspace(1, numel( LMT2AP), numel(xcom_z)), 'spline')';
        LLMAP = interp1(1:numel(LLMAP), LLMAP, linspace(1, numel( LLMAP), numel(xcom_z)), 'spline')';
        RMT2AP = interp1(1:numel(RMT2AP),  RMT2AP, linspace(1, numel( RMT2AP), numel(xcom_z)), 'spline')';
        RLMAP = interp1(1:numel(RLMAP),  RLMAP, linspace(1, numel(RLMAP), numel(xcom_z)), 'spline')';
        LMT2ML = interp1(1:numel(LMT2ML), LMT2ML, linspace(1, numel( LMT2ML), numel(xcom_x)), 'spline')';
        LLMML = interp1(1:numel(LLMML), LLMML, linspace(1, numel( LLMML), numel(xcom_x)), 'spline')';
        RMT2ML = interp1(1:numel(RMT2ML),  RMT2ML, linspace(1, numel( RMT2ML), numel(xcom_x)), 'spline')';
        RLMML= interp1(1:numel(RLMML),  RLMML, linspace(1, numel(RLMML), numel(xcom_x)), 'spline')';
    end
    %
    % x=1;% x: Left Medial (M) (1)
    % y=2;% y: Vertical Up (2)
    % z=3;% z: Forward Anterior (3)
    
    if TRCs(i).data.LMT2(i,3) > TRCs(i).data.RMT2(i,3) % MoS
        mos_AP(:,i) =  abs(LMT2AP - xcom_z);
        mos_ML(:,i) =  abs(LLMML - xcom_x);
        bos_AP(:,i) =  LMT2AP - RLMAP; %check this
        bos_ML(:,i) =  LLMML - RLMML;
    else
        mos_AP(:,i) =  abs(RMT2AP - xcom_z);
        mos_ML(:,i) =  abs(RLMML - xcom_x);
        bos_AP(:,i) =  RLMAP - LMT2AP;
        bos_ML(:,i) =  LLMML - RLMML;
    end
    % end
    names{i} = name;
    %     gaitEvents(i).name = TRCs(i).name;
    %     gaitEvents(i).data ={HSLocL, TOLocL,HSLocR, TOLocR};
    %     clear tHS HSLoc tTO TOLoc
    %clear comX com com_vel  
end

%% plotting
try
    subplot(2,1,1);hold on;plot(mos_AP);legend('mos AP')
    subplot(2,1,2);hold on;plot(mos_ML);legend('mos ML')
    
    figure;
catch
    subplot(2,1,1);hold on;legend('mos AP')
    subplot(2,1,2);hold on;legend('mos ML')
end

try
    subplot(2,1,1);hold on;plot(bos_AP);legend('bos AP')
    subplot(2,1,2);hold on;plot(bos_ML);legend('bos ML')
    title('BoS')
    
    figure;
catch
    subplot(2,1,1);hold on;legend('bos AP')
    subplot(2,1,2);hold on;legend('bos ML')
    title('BoS')
    figure;
end

try
    subplot(2,1,1);hold on;plot(COM(:,3,i));legend('COM AP')
    subplot(2,1,2);hold on;plot(XCOM(:,3,i));legend('XCOM AP')
    title('COM and XCOM')
catch
    subplot(2,1,1);hold on;legend('COM AP')
    subplot(2,1,2);hold on;legend('XCOM AP')
    title('COM and XCOM')
end

try
    out=regexp(source,'/','split');
    gaitUtils.subjName = out(end-1);
    gaitUtils.trialNames = names;
    gaitUtils.mos_AP = mos_AP;
    gaitUtils.mos_ML = mos_ML;
    gaitUtils.bos_AP = bos_AP;
    gaitUtils.bos_ML = bos_ML;
    gaitUtils.COM_VEL = COM_VEL;
    gaitUtils.XCOM = XCOM;
    gaitUtils.COM = COM;
    % [mos_AP mos_ML bos_AP bos_ML, com_vel XCOM COM]
    if saving == 1
        save('gaitUtils','gaitUtils')
    end
catch
    fprintf("failed to save the gaitUtils since some of the variables does not exist.")
end
% save('gaitEvents','gaitEvents')

%% some playign
% pp =mean(mos_AP);pp =mean(mos_ML)
% g={'W1+L+ACC.trc';'W1+L+DEC.trc';'W1+R+ACC.trc';'W1+R+DEC.trc';'W1+T+L+ACC.trc';'W1+T+L+DEC.trc';'W1+T+R+ACC.trc';'W1+T+R+DEC.trc'}
% c= categorical(g)

% bar(c,pp(6:13)) % first speed
% bar(c,pp(22:29))%seecond
% bar(c,pp(37:44))%third
% p3d = [pp(6:13);pp(22:29);pp(37:44)]
% bar3(p3d)
