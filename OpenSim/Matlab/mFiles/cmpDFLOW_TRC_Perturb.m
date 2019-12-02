%% Predict from TRCs
close all
clear 
clc
addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
trcPath = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\CB1\trcResults';cd(trcPath)
% trcPath = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\AR\DFLOW';
% This or the onebelow
% [framePerturb, BeltVelPredic,FootVelL, FootVelR] = dflowReproduceTRC_alignSignal(trcPath, -1, -1);
load BeltVelPredict.mat
% 

source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\CB1\';
[missginFiles, BeltVel, Perturb] = dflow(source, -1,-1);
close all
% Get what you can from actual DFLow data
trials = {'W1+L+ACC';'W1+L+DEC';'W1+R+ACC';'W1+R+DEC';'W1+T+L+ACC';'W1+T+L+DEC';'W1+T+R+ACC';'W1+T+R+DEC';...
    'W2+L+ACC';'W2+L+DEC';'W2+R+ACC';'W2+R+DEC';'W2+T+L+ACC';'W2+T+L+DEC';'W2+T+R+ACC';'W2+T+R+DEC';...
   'W3+L+ACC';'W3+L+DEC';'W3+R+ACC';'W3+R+DEC';'W3+T+L+ACC';'W3+T+L+DEC';'W3+T+R+ACC';'W3+T+R+DEC' };
%%
for i=1:length(trials)
    
index1 = find(strcmp({BeltVel.name}, trials(i))==1);
index2 = find(strcmp({BeltVelPredict.name}, trials(i))==1);

if  ~isempty(index1) && ~isempty(index2)
figure
hold on;plot(BeltVel(index1).data);plot(BeltVelPredict(index2).data);legend('Actual BeltSpeed','Predicted');title(trials(i))
end
end