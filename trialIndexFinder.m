function [ind, missTrials]=trialIndexFinder(trials1,trials2)

%*******************************
%Author(s): Hossein Mokhtarzadeh
%Date: 2 Sep 2019
% Decribe: How to get index of one cell of strings in another one?
%*******************************

% so length(trials1)<length(trials2)
% index1 is the index of traisl1 in trials2 i.e. 
% [ind ]=trialIndexFinder(trials1,trials2)
% trials1{1} = trials2{ind(1)}

% trials = {'W1+L+ACC';'W1+L+DEC';'W1+R+ACC';'W1+R+DEC';'W1+T+L+ACC';'W1+T+L+DEC';'W1+T+R+ACC';'W1+T+R+DEC';...
%     'W2+L+ACC';'W2+L+DEC';'W2+R+ACC';'W2+R+DEC';'W2+T+L+ACC';'W2+T+L+DEC';'W2+T+R+ACC';'W2+T+R+DEC';...
%    'W3+L+ACC';'W3+L+DEC';'W3+R+ACC';'W3+R+DEC';'W3+T+L+ACC';'W3+T+L+DEC';'W3+T+R+ACC';'W3+T+R+DEC' };

%% Example
% if we want to find the index of elemtns of one set of names in another variable
%1.  trials1 = {'W1+L+ACC';'W1+L+DEC'}
% 2. import dflow from a relevant folder
% for i=1:length(dflow.Perturb)
%  trials2{i}= dflow.Perturb(i).name
% end

%% findeing index 2
k=1;p=1;
ind =[];
missTrials=[];
for i=1:length(trials1)
    
% index1 = find(strcmp({BeltVel.name}, trials(i))==1);
% index2 = find(strcmp({BeltVelPredict.name}, trials(i))==1);
if  ~isempty(find(strcmp(trials2, trials1(i))==1))
ind(p) = find(strcmp(trials2, trials1(i))==1);
p=p+1;
else 
    missTrials(k) = trials1(i);
    k=k+1;
end
% index2 = find(strcmp({BeltVelPredict.name}, trials(i))==1);
% if  ~isempty(index1) && ~isempty(index2)
% figure
% hold on;plot(BeltVel(index1).data);plot(BeltVelPredict(index2).data);legend('Actual BeltSpeed','Predicted');title(trials(i))
% end
end


