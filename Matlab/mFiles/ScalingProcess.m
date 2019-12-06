%% for scaling do the following there is manual work here and check in OpenSim GUI :

%% First Copying relevenat files and change Mass and Height of the person in Scale.xml
% In the scalingSetup folder copy the following files (overall 8 files )
%
% gait2392_simbodyHBM.osim              markersHBM.xml
% Calibration.trc                       MU2392_Arms_Scale_MarkerSet.xml
% subject01_Scale_ScaleSet.xml          MU2392_Arms_Scale_MeasurementSet.xml
% MU2392_Arms_Scale_Tasks.xml           Scale.xml

%% TODOs
% Add mass, height and whether any coordinate needs to be fixed

function ScalingProcess(source, destination, changes)

%% Example
% source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\OpenSim_requirments\ScalingXMLs_TargetSearch\';
% destination = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\JS\';
% ScalingProcess(source, destination, changes) % changes.lock =1 means the
% foot joints is locked



%% Now Copy the relvenat Calibration.trc file for this person.

trcsource = strcat(destination,'/trcResults/');
destinationScale = strcat(destination,'/ScalingSetup/');
copyfile(source, destinationScale)
calibTRC = strcat(trcsource,'Calibration.trc');
if exist(calibTRC)
    copyfile(calibTRC, destinationScale)
else
    disp('The scaling cannot be done, there is no Calibration file in this folder')
end
%% Second go to the right Folder i.e. where ScalingSetup and the files above are:

cd(destinationScale)

%% make changes as to height, weight, joint changes, any naming etc.

%% Third run the following make sure you are in the right directory and right person being scaled
%  THIS EXECUTABLE IS DEPRECATED AND WILL BE REMOVED IN A FUTURE RELEASE.
%
%     Use opensim-cmd instead, which can do everything that this executable can.
%
%       scale -S SetupFileName -> opensim-cmd run-tool SetupFileName
%       scale -PS              -> opensim-cmd print-xml scale
%  !scale -S Scale.xml
% Pull in the modeling classes straight from the OpenSim distribution
import org.opensim.modeling.*
% scaleTool.setSubjectHeight(50)
scaleTool = ScaleTool('Scale.xml');

% scaleTool.setSubjectAge(changes.age)
% scaleTool.setSubjectHeight(changes.height)
scaleTool.setSubjectMass(changes.mass)
% change model
if changes.lock == 1
    scaleTool.getGenericModelMaker.setModelFileName('gait2392_simbodyHBMLocked.osim');
    
else
    scaleTool.getGenericModelMaker.setModelFileName('gait2392_simbodyHBM.osim')
end

geoPath = '/Applications/OpenSim 4.0/OpenSim 4.0.app/Contents/Resources/OpenSim/Geometry';
addpath geoPath;
copyfile(geoPath, destinationScale)
scaleTool.run()
% !opensim-cmd run-tool Scale.xml

%%  copy the scaled model into destination+ScaledModel folder
copyfile('subject01_simbody.osim', strcat(destination,'/ScaledModel'))

%% Fourth, drage and drop the scaled model in OpenSim i.e. subject01_scaledOnly.osim
% then So ask yourself: Does this look like a proper scaling, you can check
% the RMSE and compare it with the image if you have from the experiments.,
% check mainly out.txt file and look for errors. There might be errors that
% you don't see via runnig. So e.g.Error detected by Simbody method Xml::readFromFile(): Failed to load the Xml file 'MU2392_Arms_Scale_Tasks' with error 'Failed to open file (line=0, col=0)'.
%   (Required condition 'loadOK' was not met.)