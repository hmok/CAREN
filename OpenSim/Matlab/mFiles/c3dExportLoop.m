% ----------------------------------------------------------------------- %
% The OpenSim API is a toolkit for musculoskeletal modeling and           %
% simulation. See http://opensim.stanford.edu and the NOTICE file         %
% for more information. OpenSim is developed at Stanford University       %
% and supported by the US National Institutes of Health (U54 GM072970,    %
% R24 HD065690) and by DARPA through the Warrior Web program.             %
%                                                                         %
% Copyright (c) 2005-2018 Stanford University and the Authors             %
% Author(s): James Dunne                                                  %
%  % edited by Hossein Mokhtarzadeh for CAREN Lab                                                                       %
% Licensed under the Apache License, Version 2.0 (the "License");         %
% you may not use this file except in compliance with the License.        %
% You may obtain a copy of the License at                                 %
% http://www.apache.org/licenses/LICENSE-2.0.                             %
%                                                                         %
% Unless required by applicable law or agreed to in writing, software     %
% distributed under the License is distributed on an "AS IS" BASIS,       %
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or         %
% implied. See the License for the specific language governing            %
% permissions and limitations under the License.                          %
% ----------------------------------------------------------------------- %
% Author:James Dunne
% Modified and some additions (see the function at the end) by Hossein Mokhtarzadeh
%% James Comment on spikes @James Dunne <jjdunne@stanford.edu>
% It is typical to see COP and Tz spikes at the transitions of impact (heel strike and toe-off during gait) 
% since very small force values are used in the denominator of the COP/Tz and this results in numerical errors. 
% So, when processing forceplate data, it is typical to zero all force and moment values below some 
% threshold? motion analysis companies do this for you in their software. If you are using the OpenSim C3D reader,
% BTK is used read the C3D file and there is no thresholding being performed.
% So I do smoothing and then make cop in AP/ML to go back within the area and to the
% past cop x,y and moments >15Nm and GRFv <40N, zero.

%% Example of using the Matlab-OpenSim class 
function c3dExportLoop(path, filename, filter)
disp("hello")
disp(filter)
disp(filter(3));
%% Load OpenSim libs
import org.opensim.modeling.*
% if filter = -1 then no filtering occurs [filter (y/n) cutOff=15hz order=4]
% ex: c3dExportLoop(path, filename, [-1 15 4]) % no filter
% ex: c3dExportLoop(path, filename, [1 10 4]) % with filter
%% Get the path to a C3D file
% subjectDir = uigetdir('testData', 'Select the folder that contains the current subject data');
% subjectDir = dirIn;
% [filename, path] = uigetfile('*.c3d');
c3dpath = fullfile(path,'C3DFiles',filename);
disp(c3dpath)
%% Construct an opensimC3D object with input c3d path
% Constructor takes full path to c3d file and an integer for forceplate
% representation (1 = COP). 
c3d = osimC3D(c3dpath,0);
   
%% Get some stats...
% Get the number of marker trajectories
nTrajectories = c3d.getNumTrajectories();
% Get the marker data rate
rMarkers = c3d.getRate_marker();
% Get the number of forces 
nForces = c3d.getNumForces();
% Get the force data rate
rForces = c3d.getRate_force();

% Get Start and end time
t0 = c3d.getStartTime();
tn = c3d.getEndTime();

%% Rotate the data 
% c3d.rotateData('x',-90)
% to make it as OpenSim directon do the following(X:AP, Y: Vertical, Z: ML)
% after trasformation.
c3d.rotateData('y',90)
c3d.rotateData('z',90)


%% Get the c3d in different forms
% Get OpenSim tables
markerTable = c3d.getTable_markers();
forceTable = c3d.getTable_forces();
% Get as Matlab Structures


        
        [markerStruct, forceStruct] = c3d.getAsStructs();


%% Filter both Marker and Force Data
% use ForceStruct and markerStruct and filter in osimC3D and then convert
% it to obj.markers and obj.forces so how?

%% make nons zero in Mot files
       
%         fnames = {};ind=[];
%         for fn = fieldnames(forceStruct)'
%             fieldcontent = forceStruct.(fn{1});
%             if isstruct(fieldcontent)
%                 subfields = findnanfields(fieldcontent);
%                 if ~isempty(subfields)
%                     fnames = [fnames; strcat(fn{1}, '.', subfields)];
%                 end
%             elseif sum(sum(isnan(fieldcontent)))>=1 %|| (contains(fn{1},'p') && )
%                 fnames = [fnames; fn{1}];
%                 fieldcontent(isnan(fieldcontent)) = 0;
%                 forceStruct = setfield(forceStruct,fn{1},fieldcontent);
%                 %            A(A==0) = NaN;
%             end
%             % this is to ignore too much changes in CoP in eed to find a
%             % better way though this does nto seem to work at all so I will
%             % fix it later 
% %             if (contains(fn{1},'p')) 
% %             ind=find(abs(fieldcontent./1000)>5);
% %             if ~isempty(ind)
% %                 fnames = [fnames; fn{1}];
% %                 fieldcontent(ind) = 0;
% %                 forceStruct = setfield(forceStruct,fn{1},fieldcontent);
% %             end
% %             end
%         
%         end

%% if filter = -1 then no filtering occurs
if filter(1) ~= -1
    
    %% Filter both Marker and Force Data
    % use ForceStruct and markerStruct and filter in osimC3D and then convert
    % it to obj.markers and obj.forces so how?
    disp("filteration Butt is happening");
    
    fieldsmarker = fieldnames(markerStruct);fieldsforce = fieldnames(forceStruct);
    
    %%%%%%%%%%%%%%%%delete nans and make them zero %%%%%%%%%%%%%%
%     B = [markerStruct, forceStruct];
    
%     markerStruct = trc;
%         fnames = {};
%         for fn = fieldnames(markerStruct)'
%             fieldcontent = markerStruct.(fn{1});
%             if isstruct(fieldcontent)
%                 subfields = findnanfields(fieldcontent);
%                 if ~isempty(subfields)
%                     fnames = [fnames; strcat(fn{1}, '.', subfields)];
%                 end
%             elseif sum(sum(isnan(fieldcontent)))>=1%isnan(fieldcontent)
%                 fnames = [fnames; fn{1}];
%                 fieldcontent(isnan(fieldcontent)) = 0;
%                 markerStruct = setfield(markerStruct,fn{1},fieldcontent);
%                 %            A(A==0) = NaN;
%             end
%         
%         end
%         
%         
%            forceStruct = Forces;
%         fnames = {};
%         for fn = fieldnames(forceStruct)'
%             fieldcontent = forceStruct.(fn{1});
%             if isstruct(fieldcontent)
%                 subfields = findnanfields(fieldcontent);
%                 if ~isempty(subfields)
%                     fnames = [fnames; strcat(fn{1}, '.', subfields)];
%                 end
%             elseif sum(sum(isnan(fieldcontent)))>=1
%                 fnames = [fnames; fn{1}];
%                 fieldcontent(isnan(fieldcontent)) = 0;
%                 forceStruct = setfield(forceStruct,fn{1},fieldcontent);
%                 %            A(A==0) = NaN;
%             end
%         
%         end
        
    %%%%%%%%%%%%%%%%end above: delete nans and make them zero %%%%%%%%%%%%%%
    
    for i = 1:length(fieldsmarker)
        disp(rMarkers);
        if i~=length(fieldsmarker) % time is the last so not smoothing it
        markerStrNew.(fieldsmarker{i}) = smooth(markerStruct.(fieldsmarker{i}), filter(2), rMarkers, filter(3));
        %     markerStruct.(fieldsforce{i}) = setfield(markerStruct,fieldsforce{i},markerbuff);
        else
            markerStrNew.(fieldsmarker{i}) = markerStruct.(fieldsmarker{i}); % time is the last so not smoothing it
        end 
    end
    
    % in order to better get Smoothing GRF, it is better to get events and
    % make evetything zero when foot is in swing maybe (check this later.)
    % Anyhow I need the gait events.
    for i = 1:length(fieldsforce)
        %     markerStruct.(fieldsmarker{i}) = smooth(markerStruct.(fieldsmarker{i}), filter(2), rMakers, filter(3));
        if i~=length(fieldsforce)
           
%                         buff = forceStruct.(fieldsforce{i});
%             buff1 = buff(:,2);
%             if ((i==1 || i==4) && (any(buff1<0)))% remove or negative forces for FY1 and FY2 directin
%                 buff1(buff1<0) = 0;
%                 buff(:,2) = buff1;
%                 forceStruct.(fieldsforce{i}) = buff;
%             end

% if contains(char(fieldsforce{i}),'p') && (abs(gradient(forceStruct.(fieldsforce{i}))) >150)
%     buff=smooth(forceStruct.(fieldsforce{i}), filter(2), rForces, filter(3));
%     buff(abs(gradient(buff))<150) = 0
    
% else
        forceStrNew.(fieldsforce{i}) = smooth(forceStruct.(fieldsforce{i}), filter(2), rForces, filter(3));
% end
        %     forceStruct = setfield(forceStruct,fieldsforce{i},forcebuff);
        
%          buff = forceStrNew.(fieldsforce{i});
%             buff1 = buff(:,2);
%             if ((i==1 || i==4) && (any(buff1<0)))% remove or negative forces for FY1 and FY2 directin for the new Forces after
% %                 smoothing
%                 buff1(buff1<0) = 0; % based on Gordon et al 2018 the 40 N threshold comes then mae all forces, moment zero
%                 
%                 buff(:,2) = buff1;
%                 forceStrNew.(fieldsforce{i}) = buff;
%             end
        
        else
            forceStrNew.(fieldsforce{i}) = forceStruct.(fieldsforce{i});
        end
    end
%    
     
     %% from Gordon DFN, et al 2018. 2018;5:1?16. next threshold
% For the next step a threshold filter was applied to the ground reaction forces and moments 
% that set all values equal to zero when the vertical force was less than 40 N. This
%      forceStrNew = 
%     left side
%     if any(forceStrNew.f1(:,2) <40 )
%      ind = find(forceStrNew.f1(:,2) < 40);
%         forceStrNew.f1(ind,:) = 0;
%         forceStrNew.m1(ind,:) = 0;
%     end
% %     right side
%     if any(forceStrNew.f2(:,2) <40) 
%      ind = find(forceStrNew.f2(:,2) < 40);
%         forceStrNew.f2(ind,:) = 0;
%         forceStrNew.m2(ind,:) = 0;
%     end
%     
%     % CoP corection for larger than 2 m in x (AP), z (.5m in ML )direction 
%     CoPx = 1;CoPz = 0.5;momenty = 15;
%      if any(abs(forceStrNew.p1(:,1))./1000 >CoPx)
%      ind = find(abs(forceStrNew.p1(:,1))./1000 > CoPx);
%         forceStrNew.p1(ind,:) = 0;
% %         forceStrNew.m1(ind,:) = 0;
%      end
%     
%      if any(abs(forceStrNew.p2(:,1))./1000 >CoPx)
%      ind = find(abs(forceStrNew.p2(:,1))./1000 > CoPx);
%         forceStrNew.p2(ind,:) = 0;
% %         forceStrNew.m1(ind,:) = 0;
%      end
%     
%      if any(abs(forceStrNew.p1(:,3))./1000 >CoPz) 
%      ind = find(abs(forceStrNew.p1(:,3))./1000 > CoPz);
%         forceStrNew.p1(ind,:) = 0;
% %         forceStrNew.m1(ind,:) = 0;
%      end
%     
%      if any(abs(forceStrNew.p2(:,3))./1000 >CoPz)
%      ind = find(abs(forceStrNew.p2(:,3))./1000 > CoPz);
%         forceStrNew.p2(ind,:) = 0;
% %         forceStrNew.m1(ind,:) = 0;
%      end
%      % this is to remove those moment >15Nm in both forceplaters
%           if any(abs(forceStrNew.m1(:,2))./1000 >momenty)
%      ind = find(abs(forceStrNew.m1(:,2))./1000 > momenty);
%         forceStrNew.m1(ind,:) = 0;
% %         forceStrNew.m1(ind,:) = 0;
%      end
%     
%      if any(abs(forceStrNew.m2(:,2))./1000 >momenty) 
%      ind = find(abs(forceStrNew.m2(:,2))./1000 > momenty);
%         forceStrNew.m2(ind,:) = 0;
% %         forceStrNew.m1(ind,:) = 0;
%      end
     
    %  [markerStruct, forceStruct] = Filter_ButtW(markerStruct,15,rMarkers, forceStruct,15,rForces);
    % markerStruct = smooth(markerStruct, filter(2), rMakers, filter(3));
    % forceStruct = smooth(forceStruct, filter(2), rForces, filter(3));
    % plot(forceStruct.time,forceStruct.p1(:,3));
    
    
    %% Convert the Structs to OpenSim Tables using osimTableFromStruct.m
    
    markerTableNew = osimTableFromStruct(markerStrNew);
    forceTableNew = osimTableFromStruct(forceStrNew);
    
    
    %% Add required meta data back to opensim Table, in particular 'DataRate'.
    
    % see osimC3D for details of how to do this and write trc and mot files
    % properly
    markerTableNew.addTableMetaDataString('DataRate', num2str(rMarkers));
    markerTableNew.addTableMetaDataString('Units', 'mm');
    
else
    
    markerTableNew = markerTable;
    forceTableNew = forceTable;
        
end

% obj.forces = forceTable;
% obj.markers = markerTable;
%% Write the marker and force data to file

% Write marker data to trc file.
% c3d.writeTRC()                       Write to dir of input c3d.
% c3d.writeTRC('Walking.trc')          Write to dir of input c3d with defined file name.
% c3d.writeTRC('C:/data/Walking.trc')  Write to defined path input path.
% c3d.writeTRC('test_data_markers.trc');
% c3d.writeTRC()
% c3d.writeTRC(markerTable)  
 
 [filepath,name,ext] = fileparts(c3dpath);
% outputPath = strcat(filepath,'\',name,'.trc');
 outputPath = strcat(path,'\','trcResults\',name,'.trc');
     
TRCFileAdapter().write( markerTableNew, outputPath)
disp(['Marker file written to ' outputPath]);

% Write force data to mot file.
% c3d.writeMOT()                       Write to dir of input c3d.
% c3d.writeMOT('Walking.mot')          Write to dir of input c3d with defined file name.
% c3d.writeMOT('C:/data/Walking.mot')  Write to defined path input path.
% 
% This function assumes point and torque data are in mm and Nmm and
% converts them to m and Nm. If your C3D is already in M and Nm,
% comment out the internal function convertMillimeters2Meters()
% c3d.writeMOT('test_data_forces.mot');
% c3d.writeMOT();
% c3d.writeMOT(forceTable);
% [filepath,name,ext] = fileparts(c3dpath);
% outputPath = strcat(filepath,'\',name,'.mot');

 outputPath = strcat(path,'\','grfResults\',name,'.mot');
        import org.opensim.modeling.*
         % Get the forces table
%          forces = obj.getTable_forces();
         forces = forceTableNew;
         % Get the column labels
         labels = forces.getColumnLabels();
         % Make a copy
         updlabels = labels; 
          
         % Labels from C3DFileAdapter are f1, p1, m1, f2,...
         % We edit them to be consistent with requirements of viewing 
         % forces in the GUI (ground_force_vx, ground_force_px,...)
         for i = 0 : labels.size() - 1
            % Get the label as a string
            label = char(labels.get(i));
            % Transform the label depending on force, point, or moment
            if ~isempty(strfind(label,'f'))
                label = strrep(label,'f', 'ground_force_');
                label = [label '_v'];
            elseif ~isempty(strfind(label,'p'))
                label = strrep(label,'p', 'ground_force_');
                label = [label '_p'];
            elseif ~isempty(strfind(label,'m'))
                label = strrep(label,'m', 'ground_moment_');
                label = [label '_m'];
            end
            % update the label name 
            updlabels.set(i,label);
         end
         
         % set the column labels
         forces.setColumnLabels(updlabels)
         
         % Flatten the Vec3 force table
         postfix = StdVectorString();
         postfix.add('x');postfix.add('y');postfix.add('z');
         forces_flat = forces.flatten(postfix);
          
         % Change the header in the file to meet Storage conditions
          if forces_flat.getTableMetaDataKeys().size() > 0
              for i = 0 : forces_flat.getTableMetaDataKeys().size() - 1
                  % Get the metakey string at index zero. Since the array gets smaller on
                  % each loop, we just need to keep taking the first one in the array.
                  metakey = char(forces_flat.getTableMetaDataKeys().get(0));
                  % Remove the key from the meta data
                  forces_flat.removeTableMetaDataKey(metakey)
              end
          end
          % Add the column and row data to the meta key
          forces_flat.addTableMetaDataString('nColumns',num2str(forces_flat.getNumColumns()+1))
          forces_flat.addTableMetaDataString('nRows',num2str(forces_flat.getNumRows()));

          % zero the vertical grf<40 force_flat1 functon below
          forces_flat1 = zeroVerGRF(forces_flat);
          % moment correction if >15Nm
%           forces_flat2 = momentCorrect(forces_flat1);
%           % CoP correction if CoPx>1m and CoPz>.5m
%           forces_flat3 = copCorrect(forces_flat2);
          % Convert mm to m
          forces_flat_m  = convertMillimeters2Meters(forces_flat1);
          


 STOFileAdapter().write(forces_flat_m, outputPath)
 disp(['Forces file written to ' outputPath]);

% path = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\DF1\';
% path2 = path;
% extForces = Storage([path2 'test_data_forces.mot']);
% table =  extForces.exportToTable()
% Array() = osimTableToStruct(table)
%     extForces.print([path2 'External Forces\NMB_ExternalForce' '.mot'], 'w');

%% conver to mm
  function table_flat = convertMillimeters2Meters(table_flat)
            % Function to convert displacement forceplate measurements made
            % in millimeters to meters. This will convert point data (mm)
            % to m and Torque data (Nmm) to Nm.
            
            nForces = table_flat.getNumColumns();
            nRows  = table_flat.getNumRows();
            labels = table_flat.getColumnLabels();
            
            for i = 0 : nForces - 1
                % Find all point and torque colomns. Force columns will
                % have _v in the label, all columns that don't have this
                % character will be point and torque columns
                if ~contains(char(labels.get(i)),'v')
                    for u = 0 : nRows - 1
                        % Get the table value
                        c = table_flat.getDependentColumnAtIndex(i).get(u);
                        % set the table value
                        table_flat.getDependentColumnAtIndex(i).set(u,c/1000);
                    end
                end    
            end
           disp('Point and Torque values convert from mm and Nmm to m and Nm, respectively')
  end

%% a function to zero vertical forces<40 to zero.
  function table_flat1 = zeroVerGRF(table_flat1)
            % Function to convert displacement forceplate measurements made
            % in millimeters to meters. This will convert point data (mm)
            % to m and Torque data (Nmm) to Nm.
            
            nForces = table_flat1.getNumColumns();
            nRows  = table_flat1.getNumRows();
            labels = table_flat1.getColumnLabels();
            
            for i = 0 : nForces - 1
                % Find all point and torque colomns. Force columns will
                % have _v in the label, all columns that don't have this
                % character will be point and torque columns
                if contains(char(labels.get(i)),'vy')
                    for u = 0 : nRows - 1
                        % Get the table value
                        c = table_flat1.getDependentColumnAtIndex(i).get(u);
                        % set the table value
                        if c<40
                        table_flat1.getDependentColumnAtIndex(i).set(u,0);
                        end
                    end
                end    
            end
           disp('Zero below certain forces for vGRF,')
  end


%% make CoP to the one before if CoPx>1m and abs(CoPy) > .5 whcih are length and width of FP.
  function table_flat2 = copCorrect(table_flat2)
            % Function to convert displacement forceplate measurements made
            % in millimeters to meters. This will convert point data (mm)
            % to m and Torque data (Nmm) to Nm.
            
            nForces = table_flat2.getNumColumns();
            nRows  = table_flat2.getNumRows();
            labels = table_flat2.getColumnLabels();
            
            for i = 0 : nForces - 1
                % Find all point and torque colomns. Force columns will
                % have _v in the label, all columns that don't have this
                % character will be point and torque columns
                if contains(char(labels.get(i)),'px')
                    for u = 0 : nRows - 1
                        % Get the table value
                        c = table_flat2.getDependentColumnAtIndex(i).get(u);
                        % set the table value
                        if abs(c/1000)>1 && u>0
                        b=table_flat2.getDependentColumnAtIndex(i).get(u-1);
                        table_flat2.getDependentColumnAtIndex(i).set(u,0);
                        end
                    end
                end 
                
                   if contains(char(labels.get(i)),'pz')
                    for u = 0 : nRows - 1
                        % Get the table value
                        c = table_flat2.getDependentColumnAtIndex(i).get(u);
                        % set the table value
                        if abs(c/1000)>.5 && u>0
                        b=table_flat2.getDependentColumnAtIndex(i).get(u-1);
                        table_flat2.getDependentColumnAtIndex(i).set(u,0);
                        end
                    end
                end 
                
            end
           disp('CoP correction in AP and ML direction if they are large')
  end

%% correct Moment free one if they are large >15 in FP make it the one before
  function table_flat3 = momentCorrect(table_flat3)
            % Function to convert displacement forceplate measurements made
            % in millimeters to meters. This will convert point data (mm)
            % to m and Torque data (Nmm) to Nm.
            
            nForces = table_flat3.getNumColumns();
            nRows  = table_flat3.getNumRows();
            labels = table_flat3.getColumnLabels();
            
            for i = 0 : nForces - 1
                % Find all point and torque colomns. Force columns will
                % have _v in the label, all columns that don't have this
                % character will be point and torque columns
                if contains(char(labels.get(i)),'my')
                    for u = 0 : nRows - 1
                        % Get the table value
                        c = table_flat3.getDependentColumnAtIndex(i).get(u);
                        % set the table value
                        if abs(c/1000)>15 && u>0
                            b=table_flat3.getDependentColumnAtIndex(i).get(u-1);
                            table_flat3.getDependentColumnAtIndex(i).set(u,0);
                        end
                    end
                end
                
                
                
            end
            disp('Moment correction in Y direction if they are large')
  end


end


%% previuous ideas

% [a1, pertLocL] = max(abs(forceStruct.p1(:,1)));
% [a2, pertLocR] = max(abs(forceStruct.p2(:,1)));
%  [a1, pertLocL]  = min(forceStruct.f1(:,2)) ;% this is to get where pertub starts the next gait event Foot strike is when we consider 
% %  as iniail cntactg of the pertubation for that particualr foot
% [a2, pertLocR]  = min(forceStruct.f2(:,2));
% if contains(filename, '+L')
% pertLoc.side = 'Left';    
% pertLoc.data = pertLocL/rForces;
% else 
%     pertLoc.side = 'Right';    
%     pertLoc.data = pertLocR/rForces;
% end
% pertLoc.name = filename;
% 
% save('pertLoc','pertLoc');clear a1 a2 pertLocL pertLocR 