%% Steps to filter data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% (i) Use osimC3D.m to read a c3d file and get the data as Matlab Structs.
% (ii) Use your own methods to filter the data in the Matlab Structs. Keep the data stored as Structs. 
% (iii) Convert the Structs to OpenSim Tables using osimTableFromStruct.m
% (iv) Add required meta data back to opensim Table, in particular 'DataRate'. 
% (v) Use the TRCFileAdapter() and STOFIleApater() to write the marker and force tables, respectively, to file. 
% 
% Notes:
% - for (iv) see this example code for adding metadata to an OpenSim Table
% - for (v) open osimC3D in your matlab editor to see how to use some of these functions. TRCFileAdapter and STOFileAdapter are used internally.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (i) Use osimC3D.m to read a c3d file and get the data as Matlab Structs.
%% Load OpenSim libs
import org.opensim.modeling.*

%% Get the path to a C3D file
[filename, path] = uigetfile('*.c3d');
c3dpath = fullfile(path,filename);

%% Construct an opensimC3D object with input c3d path
% Constructor takes full path to c3d file and an integer for forceplate
% representation (1 = COP). 
c3d = osimC3D(c3dpath,1);
   
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
c3d.rotateData('x',-90)

%% Get the c3d in different forms
% Get OpenSim tables
markerTable = c3d.getTable_markers();
forceTable = c3d.getTable_forces();
nForces = forceTable.getNumColumns();
nRows  = forceTable.getNumRows();
labels = forceTable.getColumnLabels();
%% Get as Matlab Structures
[markerStruct forceStruct] = c3d.getAsStructs();
%Copy of Marker and FOrceStruct
markers = markerStruct;forces = forceStruct;
plot(forces.time,forces.p1(:,3));hold on
%% Filter both Marker and Force Data
% use ForceStruct and markerStruct and filter in osimC3D and then convert
% it to obj.markers and obj.forces so how?
display("filteration Butt is happening")
%  [markerStruct, forceStruct] = Filter_ButtW(markerStruct,15,rMarkers, forceStruct,15,rForces);
markerStruct = smooth(markerStruct, cutOff, sampleFreq, order); 
plot(forceStruct.time,forceStruct.p1(:,3));

 
 %% Convert the Structs to OpenSim Tables using osimTableFromStruct.m
 
 markerTable = osimTableFromStruct(markerStruct);
 forceTable = osimTableFromStruct(forceStruct);
%% Add required meta data back to opensim Table, in particular 'DataRate'. 

% see osimC3D for details of how to do this and write trc and mot files
% properly
markerTable.addTableMetaDataString('DataRate', num2str(rMarkers));
markerTable.addTableMetaDataString('Units', 'mm');
% forceTable.addTableMetaDataString('DataRate', num2str(rForces));
 % Add the column and row data to the meta key
 
 
 labels = forceTable.getColumnLabels();
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
         forceTable.setColumnLabels(updlabels)
         
         % Flatten the Vec3 force table
         postfix = StdVectorString();
         postfix.add('x');postfix.add('y');postfix.add('z');
         forceTable = forceTable.flatten(postfix);
         
 
 if forceTable.getTableMetaDataKeys().size() > 0
     for i = 0 : forceTable.getTableMetaDataKeys().size() - 1
         % Get the metakey string at index zero. Since the array gets smaller on
         % each loop, we just need to keep taking the first one in the array.
         metakey = char(forceTable.getTableMetaDataKeys().get(0));
         % Remove the key from the meta data
         forceTable.removeTableMetaDataKey(metakey)
     end
 end
 % Add the column and row data to the meta key 
 
 forceTable.addTableMetaDataString('nColumns',num2str(forceTable.getNumColumns()+1))
 forceTable.addTableMetaDataString('nRows',num2str(forceTable.getNumRows()));

%% Use the TRCFileAdapter() and STOFIleApater() to write the marker and force tables, respectively, to file. 
%% Write the marker and force data to file

% Write marker data to trc file.
% c3d.writeTRC()                       Write to dir of input c3d.
% c3d.writeTRC('Walking.trc')          Write to dir of input c3d with defined file name.
% c3d.writeTRC('C:/data/Walking.trc')  Write to defined path input path.
% c3d.writeTRC('test_data_markers.trc');

  % Convert mm to m
%             outputPath = generateOutputPath(markerTable,path,'.trc');
    outputPath = [path 'test_data_markers.trc'];
           
            % Write to file
            import org.opensim.modeling.*
            TRCFileAdapter().write( markerTable, outputPath)
%              TRCFileAdapter().write( obj.markers, outputPath)
            disp(['Marker file written to ' 'test_data_markers.trc']);

% Write force data to mot file.
% c3d.writeMOT()                       Write to dir of input c3d.
% c3d.writeMOT('Walking.mot')          Write to dir of input c3d with defined file name.
% c3d.writeMOT('C:/data/Walking.mot')  Write to defined path input path.
% 
% This function assumes point and torque data are in mm and Nmm and
% converts them to m and Nm. If your C3D is already in M and Nm,
% comment out the internal function convertMillimeters2Meters()
% c3d.writeMOT('test_data_forces.mot');
          % Write to file
          outputPath = [path 'test_data_forces.mot'];
          STOFileAdapter().write(forceTable, outputPath)
%           STOFileAdapter().write(forces_flat_m, outputPath)
%           STOFileAdapter().write(forces_flat_m, outputPath)
%           disp(['Forces file written to ' outputPath]);
%           disp(['Forces file written to ' outputPath]);
%%
     function outputPath = generateOutputPath(obj,path, ext)
            % Function to generate an output path from no, partial, or fully
            % defined user path. 
            
            % Validate the output filename
            if size(path,2) > 1
                % Path object should be of size == 1, any larger and user
                % input multiple variables into function. 
                error([ num2str(size(path,2)) ' inputs, expecting zero or one'])
            end
        
            if isempty(path)
               % No file path has been input, so use the path and name from
               % the c3d file. 
               filepath = obj.getPath();
               name = obj.getName();
            else
            
                if ~ischar(path{1})
                   error('Input must be a sting of characters')
                end

                if isempty(strfind(path{1}, ext))
                   error(['Input must be a path to a ' ext ' file']);
                end
            
                % User has included a path to write to
                [filepath, name, e] = fileparts(path{1}); 
                if isempty(filepath)
                  % Only the file name is given
                  filepath = obj.getPath();
                end
            end
            % Generate the output path.
            outputPath = fullfile(filepath, [name ext]);
        end