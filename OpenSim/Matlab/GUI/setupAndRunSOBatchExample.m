% ----------------------------------------------------------------------- %
% The OpenSim API is a toolkit for musculoskeletal modeling and           %
% simulation. See http://opensim.stanford.edu and the NOTICE file         %
% for more information. OpenSim is developed at Stanford University       %
% and supported by the US National Institutes of Health (U54 GM072970,    %
% R24 HD065690) and by DARPA through the Warrior Web program.             %
%                                                                         %   
% Copyright (c) 2005-2017 Stanford University and the Authors             %
%                                                                         %   
% Licensed under the Apache License, Version 2.0 (the "License");         %
% you may not use this file exceptdirectit in compliance with the License.        %
% You may obtain a copy of the License at                                 %
% http://www.apache.org/licenses/LICENSE-2.0.                             %
%                                                                         % 
% Unless required by applicable law or agreed to in writing, software     %
% distributed under the License is distributed on an "AS IS" BASIS,       %
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or         %
% implied. See the License for the specific language governing            %
% permissions and limitations under the License.                          %
% ----------------------------------------------------------------------- %
% edited by Hossein Mokhtarzadeh for CAREN Lab
% Pull in the modeling classes straight from the OpenSim distribution

function [errorCollect]=setupAndRunSOBatchExample(source, destination, rraScaled)

%% Example
% addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
% source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\OpenSim_requirments\SOXMLs_TargetSerach';
% destination = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';
% [errorCollect]=setupAndRunSOBatchExample(source, destination, 1)


%% Now Copy the relvenat xml file for this person.

ik_results_folder = fullfile(destination,'IKResults\');
GRF_results_folder = fullfile(destination,'grfResults\');
results_folder = fullfile(destination,'SOResults\');
setupfiles_folder = fullfile(destination,'SOSetup\');
% genericSetupForIK = 'IK_HBM_Setup.xml';
copyfile(source, setupfiles_folder)
modelFile = 'subject01_simbody.osim';
model_folder = fullfile(destination,'ScaledModel\');


import org.opensim.modeling.*

% move to directory where this subject's files are kept
% subjectDir = uigetdir('./DF1', 'Select the folder that contains the current subject data');

% Go to the folder in the subject's folder where IK Results are
% ik_results_folder = fullfile(subjectDir, 'IKResults');

% Go to the folder in the subject's folder where GRF data Results are
% GRF_results_folder = fullfile(subjectDir, 'grfResults');

% specify where setup files will be printed.
% setupfiles_folder = fullfile(subjectDir, 'SOSetup');

% specify where results will be printed.
% results_folder = fullfile(subjectDir, 'SOResults');

% specify where scaledModel will be printed.
% model_folder = fullfile(subjectDir, 'ScaledModel');
% %% To send an email at the end of this function define these variables appropriately:
% % more details available here:
% % http://www.mathworks.com/support/solutions/en/data/1-3PRRDV/
% %
% mail = 'youremailaddress@gmail.com'; %Your GMail email address
% password = 'yourPassword'; %Your GMail password
% 
% % Then this code will set up the preferences properly:
% setpref('Internet','E_mail',mail);
% setpref('Internet','SMTP_Server','smtp.gmail.com');
% setpref('Internet','SMTP_Username',mail);
% setpref('Internet','SMTP_Password',password);
% props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.auth','true');
% props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
% props.setProperty('mail.smtp.socketFactory.port','465');

%% Get and operate on the files
import org.opensim.modeling.*

% subjectNumber = 'subject01';
% genericSetupForAn = fullfile(setupfiles_folder, [subjectNumber '_Setup_Analyze_generic.xml']);
genericSetupForAn = fullfile(setupfiles_folder, ['\SO_HBM_Setup.xml']);

% get the file names that match the ik_reults convention
% this is where consistent naming conventions pay off
trialsForAn = dir(fullfile(ik_results_folder, '*_ik.mot'));
trialsForAnGRF = dir(fullfile(GRF_results_folder, '*.mot'));

nTrials =length(trialsForAn);
statop = AnalyzeTool(genericSetupForAn);
import org.opensim.modeling.*
tic
errorCollect = {};
for trial= 2:nTrials
  
    % get the name of the file for this trial
    motIKCoordsFile = trialsForAn(trial).name;
    motGRFCoordsFile = trialsForAnGRF(trial).name;
    
    % create name of trial from .trc file name
    name = regexprep(motIKCoordsFile,'_ik.mot','');
    nameGRF = regexprep(motGRFCoordsFile,'.mot','');
      if rraScaled ~= 1
        osimScaled = dir(fullfile(model_folder, 'subject01_simbody.osim'));
        
    else
%         name = regexprep(motIKCoordsFile,'_ik.mot','');
        osimScaled = dir(fullfile(model_folder,['rraScaledModel_' name '.osim']));
    end
    
    % get .mot data to determine time range change this to proepr API code
    % instead of Storage!
    motCoordsData = Storage(fullfile(ik_results_folder, motIKCoordsFile));
    motGRFCoordsData = Storage(fullfile(GRF_results_folder, motGRFCoordsFile));
          
% motCoords = STOFileAdapter.read(fullfile(ik_results_folder, motIKCoordsFile));
% % motCoordsData = osimTableToStruct(motCoords);
% motGRFCoords = STOFileAdapter.read(fullfile(GRF_results_folder, motGRFCoordsFile));
% motGRFCoordsData = osimTableToStruct(motGRFCoords);
    
    
    
    initial_time = motCoordsData.getFirstTime();
    final_time = motCoordsData.getLastTime();

    statop.setResultsDir(results_folder);
    model = Model(fullfile(model_folder, osimScaled.name));
    statop.setModel(model);
%      statop.setModelFileName(fullfile(model_folder, osimScaled.name));
        model_new = [model_folder modelFile];
      statop.setModelFilename(model_new)
     
%     state = model.initSystem();
% %     statop.setStatesFromMotion(state,motCoordsData,true);
     
 
    statop.setCoordinatesFileName(fullfile(ik_results_folder, motIKCoordsFile));
    statop.setStartTime(initial_time)
    statop.setFinalTime(final_time) %Needs updating
    
%     
%     ExtLoads = ExternalLoads([model_folder '\GRF.xml'],1);
    ExtLoads = ExternalLoads([setupfiles_folder 'GRF.xml'],1); % setupfiles_folder
    extloadpath = fullfile(GRF_results_folder, motGRFCoordsFile);
    ExtLoads.setDataFileName(extloadpath);
    ExtLoads.print([setupfiles_folder 'GRF1.xml'])
%     https://github.com/opensim-org/opensim-core/issues/2076 
    statop.setExternalLoadsFileName([setupfiles_folder 'GRF1.xml'])
%     idanalyzeTool.setOutputGenForceFileName([name '_ID.sto'])
%     statop.setLowpassCutoffFrequency(2)
    outfile = ['Setup_Analyze_' name '.xml'];
    statop.print(fullfile(setupfiles_folder, outfile));
%     outfileXML = fullfile(setupfiles_folder, outfile);
% Tjere is sth wrong with SO so need to cahnge teh simple !opensm works but the run here not properly     
% statop.run();
    cd(setupfiles_folder);
    eval(['!opensim-cmd run-tool ', outfile]);
    fprintf(['Performing SO on cycle # ' num2str(trial) '\n']);
    cd(results_folder);
    % rename the out.log so that it doesn't get overwritten
     copyfile(fullfile(setupfiles_folder,'out.log'), fullfile(results_folder, [name '_out.log']));
     if exist(fullfile(results_folder,'subject01v2_StaticOptimization_controls.xml'))
     copyfile(fullfile(results_folder,'subject01v2_StaticOptimization_controls.xml'), fullfile(results_folder, [name '_StaticOptimization_controls.xml']));
     copyfile(fullfile(results_folder,'subject01v2_StaticOptimization_force.sto'), fullfile(results_folder, [name '_StaticOptimization_force.sto']));
     copyfile(fullfile(results_folder,'subject01v2_StaticOptimization_activation.sto'), fullfile(results_folder, [name '_StaticOptimization_activation.sto']));
     delete(fullfile(results_folder,'subject01v2_StaticOptimization_controls.xml'));
     delete(fullfile(results_folder,'subject01v2_StaticOptimization_force.sto'));
     delete(fullfile(results_folder,'subject01v2_StaticOptimization_activation.sto'));
     errorCollect{trial} = 'Done';
     else
         errorCollect{trial} = name;
     end
%      model.delete();
%     clear model;
%      statop.delete();
%      clear statop;
%     ExtLoads.delete()
%     clear ExtLoads;
end
%sendmail(mail,subjectNumber, ['Hello! Analysis for ' subjectNumber '  is complete']);
save('errorCollect','errorCollect')
toc