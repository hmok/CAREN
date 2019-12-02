%% this is to correct the files I have measured from beltots speed etc to the right format that has 0 and 1 for all 2 min
% select 4 files txt that are relevant
% 
%% TODOs
%  1. make it in a way that all folders can be used and less manual 
% 2. make sure we have 4 or up to 8 for each trial and then define the
% folder and run this and dflow.m ...

% 000000
% 000000
% file 1
% 000000
% 000000
% 000000
% 000000
% file 2
% 000000  ------------>            file_corrected.txt to read in dflow.m
% 000000
% 000000
% 000000
% file 3
% 000000
% 000000
% file 4
% 000000
% 000000

% Author: Hossein Mokhtarzadeh, July 2019
cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB4\DFLOW\')

% [path] = uigetdir('Select the folder that contain DFLOW file .txt','*.txt');
% cd(path)
files = {'*S+T*','*W101*','*W1+P01*','*W1+T01*','*W1+T+P01*','*W201*','*W2+P01*','*W2+T01*','*W2+T+P01*',...
    '*W301*','*W3+P01*','*W3+T01*','*W3+T+P01*'};

% files = {'*S+T*','*W101*','*W1+P01*','*W1+T01*','*W1+T+P01*','*W201*','*W2+P01*','*W2+T01*','*W2+T+P01*',...
%     '*W301*','*W3+P01*','*W3+T01*','*W3+T+P01*'};
% files = {'*TW*'}
%%
for p = 5%1:length(files)
    
ExpFolderName = dir(files{p});  % '*W2+T+P01*.txt' has a proble check the timing etc
% ExpFolderName = dir('*S+T*.txt');  % '*W2+T+P01*.txt' has a proble check the timing etc
% ExpFolderName = dir('*W1+P01*.txt');  % '*W2+T+P01*.txt' has a proble check the timing etc
n=length(ExpFolderName);
fileNames = cell(n,1);
%%
for i=1:n
    fileNames{i} = convertCharsToStrings(ExpFolderName(i,1).name);
end

disp('the number of files that merged :')
files{p}
length(fileNames)
%%
[~, name, ext] = fileparts(fileNames{1});

add = 10; % 10 row of zero vicon trigger to add beginign and end of each file
for i=1:length(fileNames) %e,g. 5 % number of files for such a trial e.g. S+T has 5 and other have 4 trials. so the aim is to save one file from all 
%      struc = importdata(fileNames{i});
     a= readtable(fileNames{i});

     sA=size(a);
     vNameL = length(a.Properties.VariableNames);
     groupData = array2table(zeros(sA(1)+2*add,vNameL));
     groupData.Properties.VariableNames = a.Properties.VariableNames;
      sData=size(groupData);
%      sA-sData
     groupData(add+1:end-add,:) = a;
     deltaT = a.Time(2)-a.Time(1);
     for p=add:-1:1
     groupData.Time(p,:) =  groupData.Time(p+1,:)-deltaT;
     end
     clear p
     for p=1:add %eend times
     groupData.Time(sData(1)-add+p,:) =  groupData.Time(sData(1)-add+p-1,:)+deltaT;
     end
%      groupData.Time(add-1,:) = a.Time(1)-deltaT;
%     b.ViconTrigger(end-add:end) = 0;
     c(i).data = groupData;
     clear groupData a 
end


%%
% d1= vertcat(c(1:5).data);
% d2= vertcat(c(1).data,c(2).data,c(3).data,c(4).data,c(5).data);

%% saving and contacating etc
p=  length(fileNames);
d= vertcat(c(1:p).data);

saveName = strcat(name,'_corrected',ext);

writetable(d,saveName,'Delimiter','\t')

end
%% then run dflow for checking
% run dflow.m

% for i = 2:2:size(c.data)
% d = vertcat(c(1).data,c(2).data); 
% end
% d = vertcat(b,a); 
% 
% saveName = strcat(name,'_corrected',ext);
% 
%      writetable(b,saveName,'Delimiter','\t')


% then a.data = add 0 to end+1
% then add b.data to it
% then add add 0 to end+1
% then add c.data to it
% then add add 0 to end+1
% then add d.data to it
% then add add 0 to end+1
% 
% now you write this new file back as BG_W1_