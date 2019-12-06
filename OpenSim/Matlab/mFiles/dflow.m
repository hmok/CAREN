%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is to get DFLOw and make it as the frequncy of C3d i..e 100Hz for now
%% clean up
%% TODOs
% 1. the one that is missing i.e. missginFiles go and reproduce it from this file
% trcPath = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\VICON\CB1\trcResults';
% [framePerturb, BeltVelPredic,FootVelL, FootVelR] = dflowReproduceTRC_alignSignal(trcPath, -1, -1)
% or I can reproduce all with the above and them compare with wat i get
% here and then if ok use the above...thinka about it1


%% what it does?
% it reads the txt file from DFLOW and then finds where the belt velocities (including pertubbation happens)
% inputs: txt file from DFLOW e.g.
% outputs:
% then do the follwing:
% TODOS
% 1. make sure what the files look like and then i have all , correct the
% namings
% 2. compare with what I get from dflowReproduceTRC.m (how?)
% 3. do we need to chaeck with TRC and foot velocity
% 4.
% 5.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author(s): Dr. Hossein Mokhtarzadeh
% Date: July 30 2019
% Initial Date: July 9 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% main scipt

function [missingFiles, BeltVel, Perturb] = dflow(source, plotting, saving)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example
% make sure all files are nanmed from dflow properly and in the dflow folder then use this 
% 1. 
% addpath('C:\Users\mhossein\Documents\OpenSim\4.0\Code\Matlab');
% 2. 
% dflowDir 
% source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\CB1\';
% source = 'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\SH\';
% 3. 
% [missingFiles, BeltVel, Perturb] = dflow(source, 1,-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% source
folder = fullfile(source,'DFLOW');%'C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\TargetSearch\CB1\DFLOW';
% plotting =-1;
close all

% import org.opensim.modeling.*
% [path] = uigetdir('Select the folder that contain DFLOW file .txt','*.txt');
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB1\DFLOW')
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Projects\Data\TargetSearch\CB1\DFLOW')
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\NT4\DFLOW')
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\AN\DFLOW')
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\AH\DFLOW')
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\NT5\DFLOW')
% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\VC\DFLOW')
% cd(path)
% CB_W1+P01
freq = 100;
% vel = max()
c3dnames = {'S+T0','W10','W1+P0','W1+T0','W1+T+P0', ...
    'W20','W2+P0','W2+T0','W2+T+P0', ...
    'W30','W3+P0','W3+T0','W3+T+P0'};


cd(folder)
ExpFolderName = dir(folder);

n=length(ExpFolderName);
% fileNames = cell(n,1);
%
for i=3:n
    fileNames{i-2} = convertCharsToStrings(ExpFolderName(i,1).name);
    [~,names,~]=fileparts(fileNames{i-2});
    fileNames{i-2} = names;
end

% sort(fileNames)
k=1;h=1;missingFiles={};
for i=1:length(c3dnames)
    for j=1:length(fileNames)
        
        p(j)= contains(fileNames{j},c3dnames{i});
        if p(j)==1
            f(k).data = k;
            f(k).name = fileNames{j};
            k=k+1;
        
        end
        
        
    end
    if sum(p) ==0
         disp('This is missing: ');c3dnames{i}
         missingFiles{h} = c3dnames{i};
         h=h+1;
    end
end

% else
%             
%          if sum(p) == 0
%             disp('This is missing: ');c3dnames{j}
%         end
%         clear p 


%% to get W,  W + T
g=1; len = 700;
v = [.556 1.111 6/3.6];
    for j=1:4%length(Lvel(1,:)) % to get W1 W2 W3
        
%          if all(Lvel(:,j) - Rvel(:,j)) == 0 && all(Lvel(:,j)) >0 && ~contains(filename, '+T+') && ~contains(filename, '+P')
            BeltVel(:,g).name = strcat('W',int2str(1),'0',int2str(j));
            BeltVel(:,g).data = v(1)*ones(len,1);
            g=g+1;
            
            BeltVel(:,g).name = strcat('W',int2str(2),'0',int2str(j));
            BeltVel(:,g).data = v(2)*ones(len,1);
            g=g+1;
            
            BeltVel(:,g).name = strcat('W',int2str(3),'0',int2str(j));
            BeltVel(:,g).data = v(3)*ones(len,1);
            g=g+1;
            
%          end
        
    end
    
    for j=1:4%length(Lvel(1,:)) % W + T  % to get W1, W2 W3 + T
        
%         if all(Lvel(:,j) - Rvel(:,j)) == 0 && all(Lvel(:,j)) >0 && contains(filename, '+T+') && ~contains(filename, '+P')
            BeltVel(:,g).name = strcat('W',int2str(1),'+T0',int2str(j));
            BeltVel(:,g).data = v(1)*ones(len,1);
            g=g+1;
            
            BeltVel(:,g).name = strcat('W',int2str(2),'+T0',int2str(j));
            BeltVel(:,g).data = v(2)*ones(len,1);
            g=g+1;
            
            BeltVel(:,g).name = strcat('W',int2str(3),'+T0',int2str(j));
            BeltVel(:,g).data = v(3)*ones(len,1);
            g=g+1;
%         end
    end
    
    

%% importdata('CB2_W1+P01.txt')

for w = 1:4%n % 3 cause the first two are ., ..
    % b = importdata('CB2_W3+P01.txt');
    %  b = importdata('Subject_BG5_W3010002_corrected.txt');
    % filename = 'CB_W1.txt';
    filename = fileNames{w};%'CB_W1.txt';
    disp(filename)
    b = importdata(strcat(filename,'.txt')); name = filename; %[filepath,name,ext] = fileparts(filename);
    a = b.data;t = a(:,1);
    % fnd the first nonzero in time
    k = find(t,1);
    time = t-a(k+1,1);
    if plotting == 1
    figure;hold on; plot(time,a(:,6),'.-');plot(time,a(:,8),'.-');hold off;legend('Left','Right')
    figure;plot(time,gradient(a(:,10)))
    end
    [~, loc] = findpeaks(abs(gradient(a(:,10))));
    
    if plotting == 1
    hold on;findpeaks(abs(gradient(a(:,10))),time);plot(time,a(:,6),'.-');plot(time,a(:,8),'.-');hold off;legend('Left','Right');hold off;
    ylim([-.6 2.3]);
    end
    [~, diffLoc] = min(diff(loc)); % the freq are differnet so get the minimum diference between locations as to find the same pattern
    % lets match them with c3d and 100Hz so fit a
    % k=1;
    % for i=1:2:8
    % Lvel(:,k) = a([loc(i):loc(i+1)],6);
    %
    %  = a(loc(i):loc(i+1),8);
    % k=k+1;
    % end
    b.textdata(6);
    [~, LbeltSpClm] = find(contains(b.colheaders, 'LBeltSpeed') ==1);
    [~, RbeltSpClm] = find(contains(b.colheaders, 'RBeltSpeed') ==1);
    k=1; % this is order of trials saved
    w
    for i=1:2:8
        L =  a([loc(i)+1:loc(i+1)],LbeltSpClm); % left speed
        R = a(loc(i)+1:loc(i+1),RbeltSpClm); % right belt speed
        x =  a([loc(i)+1:loc(i+1)],1)-a([loc(i)+1],1);
        xx = 0:1/freq:round(x(end))-1/freq;
        LL = spline(x,L,xx);
        RR = spline(x,R,xx);
        Lvel(:,k) = LL';
        Rvel(:,k) = RR';
        clear RR LL xx L R
        k=k+1;
        % plot(xx,pp)
    end
    ML = mode(Lvel); % most frequent of speed
    MR = mode(Rvel);
    
    v = [.556 1.111 6/3.6];
    [W I]=min(abs(ML(1)-v)); % I become the velocty i.e. W1 or W2 etc.
    
    
     for j=1:length(Lvel(1,:)) %W +P
        
         if  contains(filename, '+P0') && ~contains(filename, '+T+')
             
             if mean(Rvel(:,j)) > v(I) +.01
                 f = '+R+ACC';
                 sp=Rvel(:,j);
                 
             elseif mean(Rvel(:,j)) < v(I) -.01
                 f = '+R+DEC';
                 sp=Rvel(:,j);
                 
             elseif mean(Lvel(:,j)) > v(I) +.01
                 f = '+L+ACC';
                 sp=Lvel(:,j);
                 
             elseif mean(Lvel(:,j)) < v(I) -.01
                 f = '+L+DEC';
                 sp=Lvel(:,j);
                 
             elseif abs(mean(Rvel(:,j)) - v(I)) <.01 && abs(mean(Lvel(:,j)) - v(I)) <.01
                 disp('cehck this out:'); filename
                 f = '+P0';
                 sp=Lvel(:,j);
                 
             end
            
            BeltVel(:,g).name = strcat('W',int2str(I),f);
            BeltVel(:,g).data = sp;
            g=g+1;
            
        end
    end
%     
    for j=1:length(Lvel(1,:)) % W + T +P
        
        if contains(filename, '+T+P0')
           if mean(Rvel(:,j)) > v(I) +.01
                 f = '+R+ACC';
                 sp=Rvel(:,j);
                 
             elseif mean(Rvel(:,j)) < v(I) -.01
                 f = '+R+DEC';
                 sp=Rvel(:,j);
                 
             elseif mean(Lvel(:,j)) > v(I) +.01
                 f = '+L+ACC';
                 sp=Lvel(:,j);
                 
             elseif mean(Lvel(:,j)) < v(I) -.01
                 f = '+L+DEC';
                 sp=Lvel(:,j);
                 
             elseif abs(mean(Rvel(:,j)) - v(I)) <.01 && abs(mean(Lvel(:,j)) - v(I)) <.01
                 disp('cehck this out:'); filename
                 f = '+P0';
                 sp=Lvel(:,j);
                 
             end
            
            BeltVel(:,g).name = strcat('W',int2str(I),'+T',f);
            BeltVel(:,g).data = sp;
            g=g+1;
            
        end
        
    end
    
   
    
      
    if plotting ==1
        figure
        subplot(2,2,1)
        hold on;plot(Lvel(:,1),'.-');plot(Rvel(:,1),'-');legend('Left','Right');tit = strcat('first in 2min, filename: ',name);title(tit);
        subplot(2,2,2)
        hold on;plot(Lvel(:,2),'.-');plot(Rvel(:,2),'-');legend('Left','Right');title('second');
        subplot(2,2,3)
        hold on;plot(Lvel(:,3),'.-');plot(Rvel(:,3),'-');legend('Left','Right');title('third');
        subplot(2,2,4)
        hold on;plot(Lvel(:,4),'.-');plot(Rvel(:,4),'-');legend('Left','Right');title('fourth');
    end
    % plot(Rvel,'--');
    % legend('Left','Right')
    % yy = ppval(pp, linspace(0,round(x(end),round(x(end)*freq)-2)));
    
end

% load framePerturb.m from  dflowReproduceTRC_alignSignal(trcPath, saving, plotting)
for jj = 1:length(BeltVel)
    [ff locc] = findpeaks((abs(gradient(BeltVel(jj).data))));
    BeltVel(jj).name
if isempty(locc)
    Perturb(:,jj).name =  BeltVel(:,jj).name;
    Perturb(:,jj).iniPerturb = 9999;
    Perturb(:,jj).endPertub = 9999;
else
    Perturb(:,jj).name =  BeltVel(:,jj).name;
    Perturb(:,jj).iniPerturb = locc(1);
    Perturb(:,jj).endPertub = locc(end);
end
end
cd ../ 
trcDir = strcat(cd,'\trcResults'); 
if exist(trcDir)
    cd(trcDir)
    if exist('predictedFrames.mat')
        load predictedFrames
        for jj = 1:length(Perturb)
        
            for ii = 1:length(predictedFrames)
                if Perturb(:,jj).iniPerturb == 9999
        
                    % fubd framePerturb the name of file and replace
                    % Perturb(:,jj).iniPerturb = framePerturb(index).data w
                    index = find(strcmp({Perturb(:,jj).name}, {predictedFrames(ii).name})==1); % PerturbFrame is the actual one frmom dflow.m
                    if ~isempty(index)
                        Perturb(:,jj).iniPerturb = predictedFrames(ii).data;
                        Perturb(:,jj).endPertub = predictedFrames(ii).data+30;
                    end
                end
        
        
            end
        end
    end
end
cd(folder)
d=2;


if plotting == 1
   disp('these figures are to check whether you collected all the files correctly or not. Chekcing once is enough!') 
    for i=1:length(BeltVel)
   figure
    plot(BeltVel(i).data)
    title(BeltVel(i).name)
end
end

if saving == 1
    dflow.missingFiles = missingFiles;
    dflow.BeltVel = BeltVel;
    dflow.Perturb = Perturb;
    save('dflow','dflow')
end
%% this is to plot all velcoties just to check so no need to invest more on this below.you can ignore

% cd('C:\Users\mhossein\OneDrive - The University of Melbourne\Summer Internships 2018\Data\CB2\DFLOW');
%
% Dflow(:,:,1) = importdata('CB2_W1+P01.txt');%first velocty
% Dflow(:,:,2) = importdata('CB2_W2+P01.txt');%2nd
% Dflow(:,:,3) = importdata('CB2_W3+P01.txt');%3rd
% v(1)=.556;v(2)=1.111;v(3)=1.667;
%
% % time = dflowW1.data(:,1)-dflowW1.data(1,1);
% LbeltSpeed = 6;RbeltSpeed = 8;
% hold on
% figure
% for i=1:3
%     time = Dflow(:,:,i).data(:,1)-Dflow(:,:,i).data(1,1);
%     Data = Dflow(:,:,i).data;
%     Max(:,:,i) = abs(max(Data));
%     subplot(3,1,i);
%     plot(time, Data(:,RbeltSpeed),'-');
%     Speed =strcat('RightBelt, ','the speed:', num2str(i));
%     title(Speed)
%     PercentR(:,:,i) = (Max(:,RbeltSpeed,i)-v(i))/v(i)*100;
% %     Percent(:,:,i) = abs(Max(:,RbeltSpeed,i)-v(i))/v(i);
% end
%
% figure
% %Left belt and foot
% for i=1:3
%     time = Dflow(:,:,i).data(:,1)-Dflow(:,:,i).data(1,1);
%     Data = Dflow(:,:,i).data;
%     Max(:,:,i) = abs(max(Data));
%     subplot(3,1,i);
%     plot(time, Data(:,LbeltSpeed),'-');
%     Speed =strcat('LeftBelt, ','the speed:', num2str(i));
%     title(Speed)
%     PercentL(:,:,i) = (Max(:,LbeltSpeed,i)-v(i))/v(i)*100;
% %     Percent(:,:,i) = abs(Max(:,RbeltSpeed,i)-v(i))/v(i);
% end


% NOw select only 1s and 0s from 10th colum Dflow.data(:,10)

% max(Dflow.data)
% plot(dflow.data(:,6),'-');
% plot(dflow.data(:,8),'.-');
% legend('LBeltS','RBeltS')
%
% LSMax=max(dflow.data(:,6));
% abs(LSMax-)
%
% plot(gradient(dflow.data(:,6)))
%     if contains(filename, 'S+T0')
%         BeltVel(:,g).name = name;
%         BeltVel(:,g).data = Rvel(:,1);
%         g=g+1;
%     elseif  contains(filename, 'W10')
%         for j=1:length(Lvel(1,:))
%             BeltVel(:,g).name = strcat('W10', int2str(j));
%             BeltVel(:,g).data = Rvel(:,j);
%             g=g+1;
%         end
%         
%     elseif  contains(filename, 'W20')
%         for j=1:length(Lvel(1,:))
%             BeltVel(:,g).name = strcat('W20', int2str(j));
%             BeltVel(:,g).data = Rvel(:,j);
%             g=g+1;
%         end
%         
%         
%     elseif  contains(filename, 'W30')
%         for j=1:length(Lvel(1,:))
%             BeltVel(:,g).name = strcat('W30', int2str(j));
%             BeltVel(:,g).data = Rvel(:,j);
%             g=g+1;
%         end
%         
%         
%     elseif  contains(filename, 'W1+T0')
%         for j=1:length(Lvel(1,:))
%             BeltVel(:,g).name = strcat('W1+T0', int2str(j));
%             BeltVel(:,g).data = Rvel(:,j);
%             g=g+1;
%         end
%         
%     elseif  contains(filename, 'W2+T0')
%         for j=1:length(Lvel(1,:))
%             BeltVel(:,g).name = strcat('W2+T0', int2str(j));
%             BeltVel(:,g).data = Rvel(:,j);
%             g=g+1;
%         end
%         
%     elseif  contains(filename, 'W3+T0')
%         for j=1:length(Lvel(1,:))
%             BeltVel(:,g).name = strcat('W3+T0', int2str(j));
%             BeltVel(:,g).data = Rvel(:,j);
%             g=g+1;
%         end
%         
%     elseif  contains(filename, 'W1+P0')
%         for j=1:length(Lvel(1,:))
%             
%             if (abs(max(Rvel(:,j))) - abs(max(Lvel(:,j)))) > .1 && max(Rvel(:,j)) > 0 % right , ACC
%                 BeltVel(:,g).data = Rvel(:,j);
%                 BeltVel(:,g).name = strcat('W1+', 'R+', 'ACC');
%                 g=g+1;
%             elseif (abs(max(Rvel(:,j))) - abs(max(Lvel(:,j)))) > .1 && min(Rvel(:,j)) < 0 % right , DEC
%                 BeltVel(:,g).name = strcat('W1+', 'R+', 'DEC');
%                 g=g+1;
%                 
%             elseif (abs(max(Lvel(:,j))) - abs(max(Rvel(:,j)))) > .1 && max(Lvel(:,j)) > 0 % left , ACC
%                 BeltVel(:,g).name = strcat('W1+', 'L+', 'ACC');
%                 g=g+1;
%                 
%             elseif (abs(max(Lvel(:,j))) - abs(max(Rvel(:,j)))) > .1 && min(Lvel(:,j)) < 0 % left , DEC
%                 BeltVel(:,g).name = strcat('W1+', 'L+', 'DEC');
%                 g=g+1;
%                 
%             else
%                 disp('this one has a proble check it:')
%                 filename
%                 
%             end
%             
%             
%         end
%             
%         
%         
%         
%     elseif  contains(filename, 'W2+P0')
%         for j=1:length(Lvel(1,:))
%             if (abs(max(Rvel(:,j))) - abs(max(Lvel(:,j)))) > .1 && max(Rvel(:,j)) > 0 % right , ACC
%                 BeltVel(:,g).data = Rvel(:,j);
%                 BeltVel(:,g).name = strcat('W2+', 'R+', 'ACC');
%                 g=g+1;
%             elseif (abs(max(Rvel(:,j))) - abs(max(Lvel(:,j)))) > .1 && min(Rvel(:,j)) < 0 % right , DEC
%                 BeltVel(:,g).name = strcat('W2+', 'R+', 'DEC');
%                 g=g+1;
%                 
%             elseif (abs(max(Lvel(:,j))) - abs(max(Rvel(:,j)))) > .1 && max(Lvel(:,j)) > 0 % left , ACC
%                 BeltVel(:,g).name = strcat('W2+', 'L+', 'ACC');
%                 g=g+1;
%                 
%             elseif (abs(max(Lvel(:,j))) - abs(max(Rvel(:,j)))) > .1 && min(Lvel(:,j)) < 0 % left , DEC
%                 BeltVel(:,g).name = strcat('W2+', 'L+', 'DEC');
%                 g=g+1;
%                 
%             else
%                 disp('this one has a proble check it:')
%                 filename
%                 
%             end
%         end
%         
%         
%     elseif  contains(filename, 'W3+P0')
%         for j=1:length(Lvel(1,:))
%             if (abs(max(Rvel(:,j))) - abs(max(Lvel(:,j)))) > .1 && max(Rvel(:,j)) > 0 % right , ACC
%                 BeltVel(:,g).data = Rvel(:,j);
%                 BeltVel(:,g).name = strcat('W3+', 'R+', 'ACC');
%                 g=g+1;
%             elseif (abs(max(Rvel(:,j))) - abs(max(Lvel(:,j)))) > .1 && min(Rvel(:,j)) < 0 % right , DEC
%                 BeltVel(:,g).name = strcat('W3+', 'R+', 'DEC');
%                 g=g+1;
%                 
%             elseif (abs(max(Lvel(:,j))) - abs(max(Rvel(:,j)))) > .1 && max(Lvel(:,j)) > 0 % left , ACC
%                 BeltVel(:,g).name = strcat('W3+', 'L+', 'ACC');
%                 g=g+1;
%                 
%             elseif (abs(max(Lvel(:,j))) - abs(max(Rvel(:,j)))) > .1 && min(Lvel(:,j)) < 0 % left , DEC
%                 BeltVel(:,g).name = strcat('W3+', 'L+', 'DEC');
%                 g=g+1;
%                 
%             else
%                 disp('this one has a proble check it:')
%                 filename
%                 
%             end
%         end
%         
%         
%     elseif  contains(filename, 'W1+T+P0')
%         for j=1:length(Lvel(1,:))
%             if (abs(max(Rvel(:,j))) - abs(max(Lvel(:,j)))) >.1 && max(Rvel(:,j)) > 0 % right , ACC
%                 BeltVel(:,g).data = Rvel(:,j);
%                 BeltVel(:,g).name = strcat('W1+','T+','R+', 'ACC');
%                 g=g+1;
%             elseif (abs(max(Rvel(:,j))) - abs(max(Lvel(:,j)))) > .1 && min(Rvel(:,j)) < 0 % right , DEC
%                 BeltVel(:,g).name = strcat('W1+','T+', 'R+', 'DEC');
%                 g=g+1;
%                 
%             elseif (abs(max(Lvel(:,j))) - abs(max(Rvel(:,j)))) > .1 && max(Lvel(:,j)) > 0 % left , ACC
%                 BeltVel(:,g).name = strcat('W1+','T+', 'L+', 'ACC');
%                 g=g+1;
%                 
%             elseif (abs(max(Lvel(:,j))) - abs(max(Rvel(:,j)))) > .1 && min(Lvel(:,j)) < 0 % left , DEC
%                 BeltVel(:,g).name = strcat('W1+','T+', 'L+', 'DEC');
%                 g=g+1;
%                 
%             else
%                 disp('this one has a proble check it:')
%                 filename
%                 
%             end
%         end
%         
%         
%         
%     elseif  contains(filename, 'W2+T+P0')
%         for j=1:length(Lvel(1,:))
%              
%             if (abs(max(Rvel(:,j))) - abs(max(Lvel(:,j)))) > .1 && max(Rvel(:,j)) > 0 % right , ACC
%                 BeltVel(:,g).data = Rvel(:,j);
%                 BeltVel(:,g).name = strcat('W2+','T+','R+', 'ACC');
%                 g=g+1;
%             elseif (abs(max(Rvel(:,j))) - abs(max(Lvel(:,j)))) > .1 && min(Rvel(:,j)) < 0 % right , DEC
%                 BeltVel(:,g).name = strcat('W2+','T+', 'R+', 'DEC');
%                 g=g+1;
%                 
%             elseif (abs(max(Lvel(:,j))) - abs(max(Rvel(:,j)))) > .1 && max(Lvel(:,j)) > 0 % left , ACC
%                 BeltVel(:,g).name = strcat('W2+','T+', 'L+', 'ACC');
%                 g=g+1;
%                 
%             elseif (abs(max(Lvel(:,j))) - abs(max(Rvel(:,j)))) > .1 && min(Lvel(:,j)) < 0 % left , DEC
%                 BeltVel(:,g).name = strcat('W2+','T+', 'L+', 'DEC');
%                 g=g+1;
%                 
%             else
%                 disp('this one has a proble check it:')
%                 filename
%                 
%             end
%         end
%         
%         
%         elseif  contains(filename, 'W3+T+P0')
%         for j=1:length(Lvel(1,:))
%           
%             if (abs(max(Rvel(:,j))) - abs(max(Lvel(:,j)))) > .1 && max(Rvel(:,j)) > 0 % right , ACC
%                 BeltVel(:,g).data = Rvel(:,j);
%                 BeltVel(:,g).name = strcat('W3+','T+','R+', 'ACC');
%                 g=g+1;
%             elseif (abs(max(Rvel(:,j))) - abs(max(Lvel(:,j)))) > .1 && min(Rvel(:,j)) < 0 % right , DEC
%                 BeltVel(:,g).name = strcat('W3+','T+', 'R+', 'DEC');
%                 g=g+1;
%                 
%             elseif (abs(max(Lvel(:,j))) - abs(max(Rvel(:,j)))) > .1 && max(Lvel(:,j)) > 0 % left , ACC
%                 BeltVel(:,g).name = strcat('W3+','T+', 'L+', 'ACC');
%                 g=g+1;
%                 
%             elseif (abs(max(Lvel(:,j))) - abs(max(Rvel(:,j)))) > .1 && min(Lvel(:,j)) < 0 % left , DEC
%                 BeltVel(:,g).name = strcat('W3+','T+', 'L+', 'DEC');
%                 g=g+1;
%                 
%             else
%                 disp('this one has a proble check it:')
%                 filename
%                 
%             end
%         end
%         
%         
%         
%     end
%     
%   % how to clean and name files?
% dflow and my prediction (trc or gait event) dont match yet! why? time
% velocity increase or decrease (is this pertubation), when foot touches
% the belt and when changes speed happnes. I think belt speed measn pertub?
% write a code to provide proper TXT files from DFLOw from the wtrong ones
% and  correct trhe dflow soon. get 4 txt files and create one siimlar to
% "CB2_W3+P01.txt"
% . clearn up all files txt read them add before and after all 4 files add
% zeros and make them one 2 min also go and check w CAREN to get a proper
% one.. Ok
% 1. name them correctly as follows: CB2_W1+P01.txt i.e.
% initials(sesseion)_c3dName.txt
% 2. find L and R and kind of Pertub and missing then save them all. then
% compare with dflowReproduceTRC.m outcomes
% 3. detect which trial R L W P T etc..
%
% make sure allt the names are correct in the DFLOW folder
% function [RBeltVel LBeltVel] = dflow(dflowDir)