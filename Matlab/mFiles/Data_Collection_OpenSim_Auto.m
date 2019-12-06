% ************************************************************
% Create Structure Array
% https://www.mathworks.com/help/matlab/matlab_prog/create-a-structure-array.html
% Author(s): Hossein Mokhtarzadeh
% Date: 27 June 2019

% Describtion: Let's get all the TRC, Mot, MoS, etc files and get one Structure file i.e. TRCs.mat
% for visualiaton later, plot,etc
% ************************************************************

% only need to select the type of format and folder to collect all relevant
% files as structure in MAT format.

% TODOs: compare Table and Struc and see how to save the data to visualize
% better later
% Make sure you collect all dta from folders so loop in folders for all subject
% can we get all the data (GRF, MOS, Muscle, etc) in one go as the folders
% are in the same place. so think about it
% ultimate goal: get all results in one MAT file which is a Structur Array.
% Do we need plotting at this stage or leave it for next step

%% Let's collect any e.g. MOT, STO, TRC files
% add funciton like below
function Data_Collection_OpenSim_Auto(source, FileType)
%source = .../targetseach/ss.;

%% Example
% %%FileType = {'trc', 'IK', 'ID', 'SO', 'RRA', 'CMC', 'FD', 'GRF',...}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FileType = {'trc', 'IK', 'ID', 'SO', 'GRF'};
% source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';
% for i=1:length(FileType)
% Data_Collection_OpenSim_Auto(source,FileType{i})
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% clear;clc

% import org.opensim.modeling.*
% addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab')
% Go to the folder in the subject's folder where .trc files are

% fileFormat = input('What File Format do you want to collect? e.g. *.trc or *.mot, etc: ','s');
%      if isempty(reply)
%         reply = 'Y';
%      end
% data_folder = uigetdir('Select the folder that contains the results you want to collect in above format.');

switch FileType %FileType = {'trc', 'IK', 'ID', 'SO', 'RRA', 'CMC', 'FD', 'GRF',...}
    case 'trc'
        data_folder = fullfile(source,'trcResults');
        fileFormat = '*.trc';
    case 'IK'
        data_folder = fullfile(source,'IKResults');
        fileFormat = '*.mot';
    case 'GRF'
        data_folder = fullfile(source,'grfResults');
        fileFormat = '*.mot';
    case 'SO'
        data_folder = fullfile(source,'SOResults');
        fileFormat = '*_force.sto';
    case 'CMC'
        data_folder = fullfile(source,'CMCResults');
        fileFormat = '*_force.sto';
    case 'ID'
        data_folder = fullfile(source,'IDResults');
        fileFormat = '*.sto'; % check this later on
    case 'RRA'
        data_folder = fullfile(source,'RRAResults');
        fileFormat = '*_force.sto';
    case 'FD'
        data_folder = fullfile(source,'FDResults');
        fileFormat = '*_force.sto'; % check this later on
end

cd(data_folder)
[filepath,name,ext] = fileparts(data_folder);
k = strfind(name,'Result'); %find Result and then go till then for the savgin
saveName=strcat(upper(name(1:k-1)),'s.mat'); % select the first two e.g. IK grf etc to save later...

% e.g. cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1')% later on ask for this or bring all the results folders for TRCs
dirIn = pwd;
% trcpath = fullfile(dirIn,'trcResults');
% cd(trcpath);
% fileFormat = ext;
ExpFolderName = dir(fileFormat);
%  ExpFolderName = dir(ext);

n=length(ExpFolderName);
fileNames = cell(n,1);
%%
for i=1:n
    fileNames{i} = convertCharsToStrings(ExpFolderName(i,1).name);
end
if strcmp(fileFormat, '*.trc')
    % lets collect all TRC files from this folder
    parfor i=1:n
        import org.opensim.modeling.*
        dirName2 = dir(fullfile(ExpFolderName(i,1).name));%dir('Test*');
        dirName2.name;
        trctimeSeriesTable = org.opensim.modeling.TRCFileAdapter.read(dirName2.name);
%         cd("/Users/duoyizhang/Intern2019Summer/OpenSim/Matlab/mFiles");
        trc = osimTableToStruct(trctimeSeriesTable);
        [filepath,name,ext] = fileparts(dirName2.name)
        TRCs(i).name = name;
        TRCs(i).data = trc;
    end
    save(saveName,'TRCs')
    clear dirIn trcpath ExpFolderName trc trctimeSeriesTable fileNames
    
    %% Let's collect Mot files or GRFs
else
    tic
    import org.opensim.modeling.*
    % cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1')% later on ask for this or bring all the results folders for TRCs
    % dirIn = pwd;
    % path = fullfile(dirIn,'grfResults');
    % cd(path);
    ExpFolderName = dir(fileFormat);
    n=length(ExpFolderName);
    fileNames = cell(n,1);
    %
    for i=1:n
        fileNames{i} = convertCharsToStrings(ExpFolderName(i,1).name);
    end
    
    % lets collect all TRC files from this folder
    for i=1:n
        
        dirName2 = dir(fullfile(ExpFolderName(i,1).name));%dir('Test*');
        dirName2.name
        if contains(dirName2.name,'Calib')
            GRFs(i).name = dirName2.name;
            continue
        end
        timeSeriesTable = org.opensim.modeling.STOFileAdapter.read(dirName2.name);
        Forces = osimTableToStruct(timeSeriesTable);
        [filepath,name,ext] = fileparts(dirName2.name);
        GRFs(i).name = name;
        GRFs(i).data = Forces;
        fields = fieldnames(GRFs(i).data);
        %     RFy = GRFs(i).data.(fields{2});
        %     time = GRFs(i).data.(fields{end});
        %     LFy = GRFs(i).data.(fields{11});
        %     gaitEve = grfGaitEvent(RFy, LFy,time, 20,-1);
        %     gaitEvent(i).name =  dirName2.name;%step length
        %     gaitEvent(i).data = gaitEve;
    end
    toc
    
    %it can be easily no switch and just use strcat in the second part of save
    %function
    if n > 0
        switch strcat(upper(saveName(1:2)))
            case 'IK'
                disp('Method is IKs')
                IKs = GRFs;
                save(saveName,'IKs')
            case 'GR'
                disp('Method is GRF')
                save(saveName,'GRFs')
            case 'SO'
                disp('Method is SOs')
                SOs=GRFs;
                save(saveName,'SOs')
            case 'RR'
                disp('Method is RRAs')
                RRAs=GRFs;
                save(saveName,'RRAs')
                
            case 'ID'
                disp('Method is IDs')
                IDs=GRFs;
                save(saveName,'IDs')
            otherwise
                disp('Unknown method.')
        end
    end
    
    % v={strcat(upper(name(1:2)))}%{'GRFs','IKs',''};
    % values= {GRFs};
    % for n=1:numel(v)
    % vn=values(n);
    % assignin('base',char(v(n)),vn{:})
    % end
    
    
    clear dirIn path ExpFolderName Forces timeSeriesTable fileNames
    
end


%% save GaitParameters
% desiredFolder = {'StepLength\','StepWidth\', 'StepTime\', 'APMoS\', 'MLMoS\'};
% gaitPara = SL;
%
%
% [LTO, LFS, RTO, RFS, timeLTO, timeLFS,timeRTO, timeRFS] = grfGaitEvent(RFy, LFy,threshold, ploting)

%ASSIGNIN Assign variable in workspace.
% v={'x1a','x2b','xA'};
% values={10,20,30};
% for n=1:numel(v)
% vn=values(n);
% assignin('base',char(v(n)),vn{:})
% end
%% For visuliazaton later on
% fields = fieldnames(GRFs(1).data)
% plot(GRFs(1).data.(fields{1}))
% GRFs(1).data


%% Lets get IKresuls from one subejct

% tic
% clear
% clc
% import org.opensim.modeling.*
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\DF1')% later on ask for this or bring all the results folders for TRCs
% dirIn = pwd;
% path = fullfile(dirIn,'IKResults');
% cd(path);
% ExpFolderName = dir('*.mot');
% n=length(ExpFolderName);
% fileNames = cell(n,1);
% %
% for i=1:2%n
%     fileNames{i} = convertCharsToStrings(ExpFolderName(i,1).name);
% end
%
% % lets collect all TRC files from this folder
% parfor i=1:n
%     dirName2 = dir(fullfile(path,ExpFolderName(i,1).name));%dir('Test*');
%     dirName2.name;
%     timeSeriesTable = org.opensim.modeling.STOFileAdapter.read(dirName2.name);
%     InKins = osimTableToStruct(timeSeriesTable);
%     IKs(i).name = dirName2.name;
%     IKs(i).data = InKins;
% end
% toc
% save('IKs.mat','IKs')
% clear dirIn path ExpFolderName Forces timeSeriesTable fileNames

