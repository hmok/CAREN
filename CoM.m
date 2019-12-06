%% Discussion in Forum
% the tests are good
% https://github.com/opensim-org/opensim-core/tree/master/Bindings/Python/tests
% Reading states in Python
% https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=91&t=9459&p=0&start=0&view=&sid=77c2d1568fcee06c57b83d54f9204db3
% load model Kinemitcs body position in Ground
% https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=91&t=8068&p=0&start=0&view=
% StatesTrajectory Class
% https://simtk.org/api_docs/opensim/api_docs/classOpenSim_1_1StatesTrajectory.html#a125b114b27974e5fa9d46039e0c63bca
% Model Class
% https://simtk.org/api_docs/opensim/api_docs/classOpenSim_1_1Model.html#a44d8bd8fb7b2aa45c118c06cfabfadb8
% API guid
% https://simtk.org/api_docs/opensim/api_docs/md_doc_APIGuide.html
% Angular Velocity and COM calcualtino my discussuion
% https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=91&t=10392&p=0&start=0&view=&sid=694a5b70e82f21d084dd030d623d801d
% common scripting commands OpenSI
% https://simtk-confluence.stanford.edu:8443/display/OpenSim/Common+Scripting+Commands#CommonScriptingCommands-BeginnerScriptingResources
%check thos out for radian to deg and vice versa 
% https://simtk-confluence.stanford.edu:8443/display/OpenSim33/Creating+Your+Own+Analysis+Part+Two

%% Matlab code to get CoM : Center of Mass using Statetrajectory 
% https://simtk.org/api_docs/opensim/api_docs/classOpenSim_1_1StatesTrajectory.html#a125b114b27974e5fa9d46039e0c63bca
%does not work yet!
clear 
clc
import org.opensim.modeling.*
model_folder = 'C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab';
model = Model(fullfile(model_folder, 'subject01_simbody.osim'));
cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\DF1\IKResults');
model.initSystem();
filename = "W101_ik.mot";
sto = Storage(filename);
% sto = STOFileAdapter.read(filename);
sto.isInDegrees() % is it in Degree
sto.setInDegrees(0) % so make it Radians I guess
% model.computeStateVariableDerivatives()
% model.setStateVariablesValue(sto)

states = StatesTrajectory.createFromStatesStorage(model, sto, true);

% timeSeriesTable = STOFileAdapter.read("W1+R+ACC_ik.mot"); 
% STOFileAdapter.write(timeSeriesTable,'test.mot');
% sto = osimTableToStruct(timeSeriesTable);
% 

% sto = model.initSystem();
% model.equilibrateMuscles(sto)
% model.equilibrateMuscles(sto);

% states = StatesTrajectory();

% model.computeStateVariableDerivatives(sto);
% model.computeStateVariableDerivatives(states)
p = model.calcMassCenterPosition(states);
% model.realizeVelocity(sto)

%% Matlab code to get CoM : Center of Mass:This one works for now.
% 
%there was a discussion here in Py: https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=91&t=10392&p=28794&start=0&view=
% this paper has good metrics stability:
% Gordon DFN, et al ., Front Robot AI. 2018;5:1–16. 
% clear 
% clc
import org.opensim.modeling.*
model_folder = 'C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab';
model = Model(fullfile(model_folder, 'subject01_simbody.osim'));
state = model.initSystem();
cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\DF1\IKResults');
motion_storage = Storage("W101_ik.mot");%MotStorage = motion_storage;
motion_storage.getSize();
n_coord = model.getNumCoordinates();

% ## Realize the state
for n_frame = 0:motion_storage.getSize()-1
      
    for i=0:n_coord-1
        
        name = model.getCoordinateSet().get(i).getName();
        index = motion_storage.getStateIndex(name);
        vector = motion_storage.getStateVector(n_frame).getData().get(index);
        model.updCoordinateSet().get(i).setValue(state,degtorad(vector));
%         speed = diff(vector);
%         model.updCoordinateSet().get(i).setSpeedValue(state,degtorad(speed))
     end
%     state = model.initSystem();
%     model.computeStateVariableDerivatives(state);
%     state = model.initSystem();
    p = model.calcMassCenterPosition(state);
%     v = model.calcMassCenterVelocity(state);
    %     COMLoc(i+1,:,n_frame+1) = str2num(char(ArrayDouble.getValuesFromVec3(p)));
    COMLoc(n_frame+1, :) = str2num(char(ArrayDouble.getValuesFromVec3(p)));
%     COMVel(n_frame+1, :) = str2num(char(ArrayDouble.getValuesFromVec3(v)));
    clear p
end

hold on; subplot(2,1,1);plot(COMLoc);title('CoMLoc'); 
% subplot(2,1,2);plot(COMVel);title('CoMVelocity');
% c = sqrt(sum((COMLoc.^2),2));
% plot(c);
% model.computeStateVariableDerivatives(state)
%
%
% model.setUseVisualizer(1);
% state = model.initSystem();
% model.getVisualizer().show(state);
%

%% center of mass  using BodyKinematics analysis
% TODO it is still incomplete so work late

import org.opensim.modeling.*

% move to directory where this subject's files are kept
subjectDir = uigetdir('testData', 'Select the folder that contains the current subject data');

% Go to the folder in the subject's folder where IK Results are
ik_results_folder = fullfile(subjectDir, 'IKResults');

% Go to the folder in the subject's folder where GRF data Results are
GRF_results_folder = fullfile(subjectDir, 'GRFResults');

% specify where setup files will be printed.
setupfiles_folder = fullfile(subjectDir, 'AnalyzeSetup');

% specify where results will be printed.
results_folder = fullfile(subjectDir, 'AnalyzeResults');

subjectNumber = 'subject01';
% genericSetupForAn = fullfile(setupfiles_folder, [subjectNumber '_Setup_Analyze_generic.xml']);
setupfiles_folder = 'C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab';
genericSetupForAn = fullfile(setupfiles_folder, ['\bodyKinematic.xml']);
analyzeTool = AnalyzeTool(genericSetupForAn);

% get the file names that match the ik_reults convention
% this is where consistent naming conventions pay off
trialsForAn = dir(fullfile(ik_results_folder, '*_ik.mot'));

% for simplicity
% trialsForAn = 
model_folder = 'C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab';
model = Model(fullfile(model_folder, 'subject01_simbody.osim'));


cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\DF1\IKResults');
% motion_storage = Storage("W301_ik.mot");
motIKCoordsFile = 'W301_ik.mot'

trialsForAnGRF = dir(fullfile(GRF_results_folder, '*.mot'));
nTrials =length(trialsForAn);

for trial= 1:nTrials
    % get the name of the file for this trial
    motIKCoordsFile = trialsForAn(trial).name;
    motGRFCoordsFile = trialsForAnGRF(trial).name;
    
    % create name of trial from .trc file name
    name = regexprep(motIKCoordsFile,'_ik.mot','');
    nameGRF = regexprep(motGRFCoordsFile,'.mot','');
    
    % get .mot data to determine time range
    motCoordsData = Storage(fullfile(ik_results_folder, motIKCoordsFile));
    motGRFCoordsData = Storage(fullfile(GRF_results_folder, motGRFCoordsFile));
    
    % for this example, column is time
        initial_time = motCoordsData.getFirstTime();
        final_time = motCoordsData.getLastTime();
    
    analyzeTool.setName(name);
    analyzeTool.setResultsDir(results_folder);
    analyzeTool.setCoordinatesFileName(fullfile(ik_results_folder, motIKCoordsFile));
    analyzeTool.setInitialTime(initial_time);
    analyzeTool.setFinalTime(final_time);   
    
    outfile = ['Setup_Analyze_' name '.xml'];
    analyzeTool.print(fullfile(setupfiles_folder, outfile));
    
    analyzeTool.run();
    fprintf(['Performing IK on cycle # ' num2str(trial) '\n']);
    
    % rename the out.log so that it doesn't get overwritten
    copyfile('out.log', fullfile(results_folder, [name '_out.log']));
    
end
%sendmail(mail,subjectNumber, ['Hello! Analysis for ' subjectNumber '  is complete']);



%% Center of mass using 4 markers


% you can easily do this:

import org.opensim.modeling.*
trctimeSeriesTable = TRCFileAdapter.read('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\DF1\trcResults\W1+R+ACC.trc');
trc = osimTableToStruct(trctimeSeriesTable);
fields = fieldnames(trc);

% Get CoM from 4 hip markes ASIS, PSIS left and right average of these markers
% LASIS= fields{4}, RASIS= fields{5},LPSIS= fields{6},  RPSIS= fields{7},
fields = fieldnames(trc);
time = trc.('time');
% Get marker data using pelvic markers
for i = 4:7
buf(:,:,i-3) = trc.(fields{i});
end
com = mean(buf,3);
figure
plot(com./1000); legend("X: ML","Y: Vertical","Z: AP")

% subplot(2,1,1);
% plot(com./1000) 
% subplot(2,2,1);plot(COMLoc);title('CoMLoc')
%
%
%% using STOFileAdapter.read
% Loading Model Kinematics with OpenSim 4.0 MATLAB API
%https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=91&t=8068&p=0&start=0&view=


% Import opensim classes
import org.opensim.modeling.*

model_folder = 'C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab';

state = model.initSystem();
cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\DF1\IKResults');
motion_storage = Storage("W1+R+ACC_ik.mot");MotStorage = motion_storage;


% Read in model from file
% model_file = 'my_file.osim';
model_file = Model(fullfile(model_folder, 'subject01_simbody.osim'));
model = Model(model_file);

% Load in the .mot file output
mot_file = IK_Results.mot';
% Copy  as .txt so can be read as table
[path, name, ext] = fileparts(mot_file);
mot_txt_file = fullfile(path, [name, '.txt']);
copyfile(mot_file, mot_txt_file);
mot_table = readtable(mot_txt_file);
delete(mot_txt_file);

% Get a reference to the state
state = model.initSystem();

% Get a reference to all bodies in model
body_set = model.getBodySet();
num_bodies = body_set.getSize;

% Get a reference to all coordinates in model
coordinate_set = model.getCoordinateSet();
num_coordinates = coordinate_set.getSize;

% Initialize positions and rotations for each body
p = cell(length(time), num_bodies);
r = cell(length(time), num_bodies);

% Loop through time
time = mot_table.time;
for time_idx = 1:length(time)
    
    % Loop through each coordinate and update position
    num_coordinates_from_mot = size(mot_table,2) - 1; % First one is time
    assert(num_coordinates_from_mot == num_coordinates, 'Check that model matches output. Number of coordinates from model does not match number of coordinates in .mot file');
   
    for coordinate_idx = 1:num_coordinates;
        
        coordinate = coordinate_set.get(coordinate_idx - 1); % idx starts with zero
        coordinate_name = char(coordinate);
        % Pull coordinate position in local frame from .mot file
        coordinate_value_local = mot_table.(coordinate_name)(time_idx);
        % Change the coordinate position in local frame
        coordinate.setValue(state, coordinate_value_local);   
    
    end
    
    % Loop through bodies and pull position in ground frame
    for body_idx = 1:num_bodies
       
        body = body_set.get(body_idx - 1); % idx starts with zero
        this_p = body.getPositionInGround(state);
        this_r = body.getTransformInGround(state).R;
        
        % Update arrays
        p{time_idx, body_idx} = this_p;
        r{time_idx, body_idx} = this_r;
    end
end


%
%
%
%
%
% manager = Manager(model);
%
% % from https://simtk.org/api_docs/opensim/api_docs/classOpenSim_1_1Manager.html#details
% state.setTime(0.0);
% manager.initialize(state);
% manager.initialize(state)
%
% dTime = 2.0;
% finalTime = 10.0;
% n = (finalTime/dTime);
% for i = 1:n
%     state = manager.setStateStorage(motion_storage);
% end
% manager.setStateStorage(motion_storage);
%
% s = model.initSystem();
% s = state
% % returns Vec3()
% p = model.calcMassCenterPosition(state);
% v = model.calcMassCenterVelocity(s);
% a = model.calcMassCenterAcceleration(s);
%
% COMLoc = str2num(char(ArrayDouble.getValuesFromVec3(p)));
% COMVel = str2num(char(ArrayDouble.getValuesFromVec3(v)));
%
%
%
%
% model.setUseVisualizer(1);
% state = model.initSystem();
% model.getVisualizer().show(state);
%
% manager = Manager(model);
% manager.setStateStorage(motion_storage);
%
% model.buildSystem()
% state = model.initializeState()
%
% motion_storage = Storage("W1+R+DEC_ik.mot")
% model.updCoordinateSet(motion_storage)
% % # Acquiring the transformation for a specific body in the bodyset already works
% model.getBodySet().get(0).getTransformInGround(state)
%
% Des_StoreData =motion_storage
%
% Des_Time = ArrayDouble();
% Des_Size = Des_StoreData.getTimeColumn (Des_Time);
%
% for i = 1:num_Coordinates
%     coordvalue = ArrayDouble();
%     Des_StoreData.getDataColumn(char(Coordinates_Name{i,1}),coordvalue);
%     coordvect = coordvalue.getAsVector;
%
%     motType = model.getCoordinateSet().get(Coordinates_Name{i,1}).getMotionType();
%
%     % Check unit for "Rotational" or "Coupled" coordinate in your case: rad or degree
%     if strcmp(motType,'Rotational')
%         Des_Data_Array(i,:) = deg2rad(osimVectorToArray(coordvect));
%     else
%         Des_Data_Array(i,:) = osimVectorToArray(coordvect);
%     end
% end
%
% t = target_time;
% for j = 1:num_Coordinates
% 	osimModel.updCoordinateSet.get(Coordinates_Name{j,1}).setValue(osimState,Des_Data_Array(j,t));
% end