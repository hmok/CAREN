%% to collect all the subject data in a simple table format in maltab

% Author: Hossein Mokhtarzadeh
% Date
% change the directory where the csv files are
%% ToDOs
% 0. trhere is a probelm with time date see what the bug is time date is
% wrong
% 1. make it automated for all namign properly SH, etc and probaly need all the names as an input
% 2. then correct HZ below this should be something generic and then saving must chage to the right subject name above
% 3. function it simialr to below wherr the source could point to excel
% file so creat one excel file with all the names as sheets then this is
% amazing...
 function [allSubjData] = Collect_Subject_Data(source)

 %% Example
% addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
 %1. change source
% source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\Publications\Biomecognition_PaperMethod\Results';
% 2. make sure this is correct eblow filename = fullfile(source,'SubjectDataFormatAll.xlsx')
% 3. [allSubjData] = Collect_Subject_Data(source)

%% conversion to varibel
% X = 'omega';
% eval([X '= [2 3 4 5]']);
%% 

filename = fullfile(source,'SubjectDataFormatAll.xlsx');
[~, description, ~] = xlsfinfo(filename);
% clear 
% clc
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\General Training\TargetSearchSubjectData')
% filename ='SubjectDataFormatAll.xlsx'; % change this file accordingly

for i=1:2%length(description)
T = readtable(filename, 'Sheet',string(description(i))); % read the file as Table in Maltab

b = char(description(i));
eval([b '= []'])% this creatte SH = [0] or the name of sheet or subject great idea.
%% Add subject data including Name, height, weight, gender, age , date and time
%e.g. Subject Name: Hongxu Zhu Height: 177 cm Weight: 70 kg Gender: Male Age: 27

HZ.name = T.nameGender(1);
HZ.height = T.height(1); %m
HZ.weight = T.weight(1); %kg
HZ.gender = T.nameGender(2); % m or f
HZ.age = T.age(1); %yrs
HZ.DNDleg = T.D_NDLeg(1); %yrs
HZ.DNDeye= T.D_NDEye(1); %yrs
HZ.gamer = T.Gamer(1); %yrs
% HZ.date = 'Feb 22 2019'; %date of test
d1 = (T.date1(1));%'2019-02-22 14:39';
d2 = (T.date2(1));%'2019-02-22 15:37';
t1 = datetime(d1,'InputFormat','yyyy-MM-dd HH:mm');
t2 = datetime(d2,'InputFormat','yyyy-MM-dd HH:mm');
HZ.iniTime = t1; %begin test
HZ.endTime = t2; %end test
HZ.note1 =T.notes(1);
HZ.note2 = T.notes(2);
HZ.data = T(:,1:11);
% SH = HZ;
eval([b '= HZ'])
save(b)
allSubjData(i).name = description(i);
allSubjData(i).data = eval(b); % THERE IS SOMETHIGN WRONG HERE



end
save('allSubjData','allSubjData')
% save('SH','SH')