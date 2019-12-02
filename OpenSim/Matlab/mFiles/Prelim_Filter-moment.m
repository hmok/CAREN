%% import
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
%    dirIn = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH'; % or loop 
  dirIn = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH';
runOpensim(dirIn,1,10,4)


%% scale
% TODOs make it a loop for all the subjects
% time = each sesseion about 2hr
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab')
% 
 source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\OpenSim_requirments\ScalingXMLs_TargetSearch\';
 destination = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';
% for i=1:n
%     %find the revelant destination and then run ScalingProcess
% end
% but I guess it is still good to do it one by oone and check all the details
changes.mass = 65;% add mass of the person kg
changes.lock = 1; % if 1 then the feet joints are locked other wise, all is free
ScalingProcess(source, destination, changes)


%% Inverse Kinematics
% time = each sesseion about ~5min for 50files or ?
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab')
 source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\OpenSim_requirments\IKXMLs_TargetSerach\';
 destination = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';

setupAndRunIKBatchExample(source, destination)

%% Inverse Dynamics
% time = each sesseion about ~5min for 50files or ?
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab')
source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\OpenSim_requirments\IDXMLs_TargetSerach\';
destination = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';

setupAndRunIDBatchExample(source, destination)

%% ploting
% grf%
grf = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\grfResults'
 g = importdata(fullfile(grf,'W1+R+ACC.mot'))
subplot(4,1,1)
plot(g.data(:,2:4),'DisplayName','g.data(:,2:4)')
title('grfs') 
%  cops
subplot(4,1,2)
plot(g.data(:,5:7),'DisplayName','g.data(:,5:7)')
title('cop')
%moments 1
subplot(4,1,3)
plot(g.data(:,8:10),'DisplayName','g.data(:,8:10)')
title('free moment')


%momoents

grf = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\IDResults'
g = importdata(fullfile(grf,'W1+R+ACC_ID.sto'))
subplot(4,1,4)
plot(g.data(:,17:18),'DisplayName','g.data(:,8:10)')
title('knee moments')


grf = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\SOResults'
g = importdata(fullfile(grf,'W1+R+ACC_StaticOptimization_activation.sto'))
subplot(4,1,4)
plot(g.data(:,17:18),'DisplayName','g.data(:,8:10)')
title('knee moments')