function [LHeelStrikes,RHeelStrikes] = HeelStrikeDetection(LFy, RFy, Cutoff, SaveName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HEELSTRIKEDETECTION:
% from https://github.com/csu-hmc/Human_Impedance_Control/blob/master/gait_average/HeelStrikeDetection.m
% Detect heel strike time points from vertical GRFs.

% Inputs: LFy (nx1 vector): Left vertical ground reaction force.
%         RFy (nx1 vector): Right vertical ground reaction force.
%         Cutoff (scalar): A cutoff value which decided whether
%         the foot touchs ground or not. 
%         SaveName (string): File name for saving generated data. 
%            if None, don't save files. 

% Outputs: LHeelStrikes (px1 vector): Left heel strike point index.
%          RHeelStrikes (qx1 vector): Right heel strike point index.

% Creator: Huawei Wang
% Time: May 21, 2019
% Location: CUS HMC lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
   case 3
      SaveName = 'None';
   case 4
   otherwise
      error('Only accept three or four inputs')
end

% Check whether LFy and RFy have the same length
if (length(LFy) ~= length(RFy))
    error('Not equal length of LFy and RFy.')
end


l = length(LFy);  % get data length
Lsign = zeros(l, 1);  % predefine sign vector as zeros
Rsign = zeros(l, 1);  % predefine sign vector as zeros

% Detect the heel strike points. Change sign vector value to 1,
% if the vertical GRF is larger than Cutoff value.
for i = 1:l
    
%     https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2384115/
% so "Heel strike was calculated at the first frame at which the vertical
% GRF was greater than 20 Newtons for a width of at least 40 samples."
% Similarly, toe off was determined at the first frame that was less than
% 20 Newtons for a width of 40 frames."
    
    if (diff(LFy) > Cutoff) && (LFy(i:i+40) > Cutoff) %
        Lsign(i) = 1;
    end
    
%     if RFy(j) > Cutoff
%         Rsign(i) = 1;
%     
%     end

end

% Generate the index of heel strike points. 
LHeelStrikes = find(diff(Lsign)==1) + 1;
RHeelStrikes = find(diff(Rsign)==1) + 1;

% Save heel strike information into files
if ~ strcmp(SaveName, 'None')
    
    fileID = fopen(strcat(SaveName, '_LHeelStrikes.txt'),'w');
    fprintf(fileID,'%d\n',round(LHeelStrikes));
    fclose(fileID);
    
    fileID = fopen(strcat(SaveName, '_RHeelStrikes.txt'),'w');
    fprintf(fileID,'%d\n',round(RHeelStrikes));
    fclose(fileID);
   
end


end
