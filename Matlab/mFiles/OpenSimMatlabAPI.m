%% Matlab OpenSim Visualizer (does not work yet!)
% TODO: does nto work yet so how can we load states (IK or even markers,
% mot) and visualize them using Matlab or Py?
import org.opensim.modeling.*
model_folder = 'C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab';
model = Model(fullfile(model_folder, 'subject01_simbody.osim'));
state = model.initSystem();
cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\DF1\IKResults');
motion_storage = Storage("W101_ik.mot");%MotStorage = motion_storage;
motion_storage.getSize();
% model.setUseVisualizer(1);
% state = model.initSystem();
% model.getVisualizer().show(state);
model.computeStateVariableDerivatives(state)

%%
