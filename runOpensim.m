%% %% run Scale, IK, ID, SO, GRF dont seem to be working in ID and SO
% I locked mtp, subtalar, and lumbar DOFs be mindful of it. give to others
% to make it more effective as we go
% edited by Hossein Mokhtarzadeh for CAREN Lab
% TODOs

%2. clean all the c3ds and dflow in proper files and folders
%3. run all together and save trc, mot, MoS, etc in the same folder as c3d
%folder sit each trial in fact
%4. then collect all the data and write in the same folder e.g. all TRCs in
%trcResults folder, etc
%5. now we hve all the results in the folders and just posprocessing is
%needed
%6. ask proper question and read what you need from the relevant folders
%and do plot, table, etc...
% 7. or we can collect all the results from all the trials in a bit STRUCT
% mat file and then just load this large data and ... ask Andres how to
% compress this..

function runOpensim(dirIn,filter,cutoff,order)

%% Example
% dirIn = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\CB1'; % or loop 
% runOpensim(dirIn,1,10,4)
%% filteration


% if filter = -1 then no filtering and nans zeroing occurs [filter (y/n) cutOff=15hz order=4]
% ex: c3dExportLoop(path, filename, [-1(No) 15 4])
% clear; clc

% function = CAREN(Dir, ..., fiklter, )
% filter = 1; 
% cutoff = 10; 
% order = 4;
% which_analyses = [1 2 3]; % 1: only trc, 2: only mot GRF, 3: both
which_analyses = 1; %to get TRCs for e.g. gait event this is easier to do it faster.
%
addpath('/Applications/OpenSim 4.0/OpenSim 4.0.app/Contents/Resources/OpenSim/Resources/Code/Matlab'); % add this path for M files
% get the folder
% dirIn = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\DF1\';
% dirIn = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB2\VICON\CB2';


% dirIn = uigetdir('Select the folder that contain dflow data and parent of c3d folder  inclduign the folder for file.c3d','*.c3d');


% dirIn = 'C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab';
dirOut = dirIn; % change if need tobe in another folder
cd(dirIn);
errCollect.name = []; %the errors in reading c3d and exporing to OpenSim



import org.opensim.modeling.*
c3dpath = fullfile(dirIn,'C3DFiles');
cd(c3dpath);
ExpFolderName = dir('*.c3d');
% eval('cd ..\..'); dirIn = cd;
cd(dirIn)
n=length(ExpFolderName);
fileNames = cell(n,1);

% Create folders for analysis
folders = {'CMCResults','RRAResults','RRASetup','AnalyzeSetup','ScalingSetup','ScaledModel','IKSetup','IKResults'...
    ,'IDSetup','IDResults','SOSetup','trcResults','SOResults','C3DFiles','grfResults','Gait Parameters'};

for i = 1:length(folders)
    if ~exist(folders{i})
        mkdir(folders{i})
    end
    
end
%%
for i=1:n
    fileNames{i} = convertCharsToStrings(ExpFolderName(i,1).name);
end
tic
%% lets do the c3d export to OpenSim
k=1;%for errCollec
% dd = [6;7;8;9;10;12;13;14;15;16;17;18;19;33;46;52;53];

for i=1:n%n%%[dd(1):dd(end)]%1:n%22:n
    % dirName2 = dir(fullfile(cd,ExpFolderName(i,1).name,'*SE*.csv'));%dir('Test*');
    dirName2 = dir(fullfile(c3dpath,ExpFolderName(i,1).name));%dir('Test*');
    dirName2.name;
    % forceFile = fullfile(cd,ExpFolderName(i,1).name,dirName2(1,1).name);%,'*.csv')
    % file1 = dir(forceFile);% file1.name
    % run('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab\c3dExport.m');
    path = dirIn;
    filename = dirName2.name;
    % cd('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab'); % add this path
    % if filter = -1 then no filtering occurs [filter (y/n) cutOff=15hz order=4]
    % ex: c3dExportLoop(path, filename, [-1 15 4]), filter = 1; cutoff = 10; order = 4;
    disp('Processing:'); %filename
    cd('/Users/duoyizhang/Intern2019Summer/OpenSim/Matlab/mFiles');
%     try
        c3dExportLoop(path, filename, [filter cutoff order])
        % EMG..import
         cd('/Users/duoyizhang/Intern2019Summer/OpenSim/Matlab/mFiles');
%      catch 
%          if warning(['Illegal Column label. *32 changed to _32'])
%            errCollect(k).name = dirName2.name;% c3dp(k).name; %Store the files names with error it also finalizes the script
%          errCollect(k).data = i;% c3dp(k).name; %Store the files names with error it also finalizes the script
%            k=k+1;
%          
%          end
        continue; 
    end
    % run('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab\c3dExportLoop.m');
    
    % !scale -S ScaleHBM_Setup.xml
    % !ik -S IK_HBM_Setupv3.xml
    % !opensim-cmd run-tool ID_HBM_Setup.xml
    % !opensim-cmd run-tool SO_HBM_Setupv3.xml
end

% save('errCollect','errCollect')


% toc