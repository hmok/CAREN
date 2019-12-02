function trials = trialsList(fileExt, fileExtPath)

% Example
% path = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\C3DFiles';
% path = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\grfResults';
% if path is given:
% trials = trialsList('*.c3d',path);
% else:
% trials = trialsList('*.c3d');

% anotehr example for GRFs
% path = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\grfResults';
% trials = trialsList('*.mot',path);

% Author: Hossein Mokhtarzadeh
% Date: July 3rd 2019

% trials to get from c3d or any other format sto, mot, trc,etc folder

%% finding the TO and HS using TRC files
% addpath('C:\OpenSim 4.0\Code\Matlab\');
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
import org.opensim.modeling.*
if nargin == 1
STR = strcat('Select the folder that contain file',fileExt);
[path] = uigetdir(STR,fileExt);
else
    [path] = fileExtPath;
end

%% let's loop for TRCs
trials  = {};
c3dp = dir(fullfile(path,fileExt));
% c3dpath = fullfile(path,c3dp(k).name);
for i=1:length(c3dp)
trials{i} = char(c3dp(i).name);
c3dp(i).name;
end
