%Read TRC and plot any markers
%% you can easily do this:

import org.opensim.modeling.*
trctimeSeriesTable = TRCFileAdapter.read('HZ5_W1+L+ACC.trc');
trc = osimTableToStruct(trctimeSeriesTable);
fields = fieldnames(trc);
figure
for i = 1:length(fields)-1
hold on
    plot(trc.(fields{i}))
end
% Left Foot 12:16, average X direction
figure
hold on;plot(trc.(fields{5}));plot(trc.(fields{6}))
legend(fields{5},fields{6})

% check if there is nan in the trc files
fnames = findnanfields(trc);

%% Get CoM from 4 hip markes ASIS, PSIS left and right average of these markers
% LASIS= fields{4}, RASIS= fields{5},LPSIS= fields{6},  RPSIS= fields{7},
fields = fieldnames(trc);
time = trc.('time');
% Get marker data using pelvic markers
for i = 4:7
buf(:,:,i-3) = trc.(fields{i});
end
com = mean(buf,3);
plot(com)

%% or do the following reading everycolumn and every elemtn that need to be done sometimes.

import org.opensim.modeling.*
model_folder = 'C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab';
% model = Model(fullfile(model_folder, 'subject01_simbody.osim'));
% state = model.initSystem();

motion_storage = Storage("HZ5_W1+L+ACC.trc");MotStorage = motion_storage;
motion_storage.getSize();
motion_storage.getColumnLabels(); % labels of colums
% n_coord = model.getNumCoordinates(); %no need if just markers or mot
% files

% ## Realize the state
for n_frame = 0:MotStorage.getSize()-1
    
    % ## initialize state
    % state = model.initSystem();
    
    for i=0:n_coord-1
        
        name = model.getCoordinateSet().get(i).getName();
        index = MotStorage.getStateIndex(name);
        vector = MotStorage.getStateVector(n_frame).getData().get(index);
        model.updCoordinateSet().get(i).setValue(state,vector);
        
%         
        %     p = model.calcMassCenterPosition(state);
        % %     COMLoc(i+1,:,n_frame+1) = str2num(char(ArrayDouble.getValuesFromVec3(p)));
        %     COMLoc(n_frame+1, :, i+1) = str2num(char(ArrayDouble.getValuesFromVec3(p)));
        
    end
    model.computeStateVariableDerivatives(state);
    p = model.calcMassCenterPosition(state);
    %     COMLoc(i+1,:,n_frame+1) = str2num(char(ArrayDouble.getValuesFromVec3(p)));
    COMLoc(n_frame+1, :) = str2num(char(ArrayDouble.getValuesFromVec3(p)));
    
end
