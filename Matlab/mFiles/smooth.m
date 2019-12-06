% smooth.m
% Smooth column data using a 4th order lowpass.
% Usage: Y = smooth(data, cutOff, sampleFreq, order) 
% Only one cut off frequency provided
% Ajay Seth
% Modified by Tim Dorn June 2009

function Y = smooth(data, Wc, sFreq, order)

%to ignore nans careful of beging and end
% the issue is the begining and last elemetn if it is Nans

% tf = isnan(data);
% ix = 1:numel(data);
%     data(tf) = interp1(ix(~tf),data(~tf),ix(tf));

% fillmissing(data,'linear');  %no good results
%   fillmissing(data, 'movmean' ,4); % no good
data = fillmissing(data,'linear','EndValues','nearest');
%   fillmissing(data,'nearest');  %no good results

if nargin == 3
    order = 4;
elseif nargin ~= 4
   error('Smooth Error: USAGE: Y = smooth(data, Wc, sFreq, order*)\n');
end

maxF = sFreq/2;

wn = Wc/maxF;

if Wc <= 0,
    Y = data;
    return
end

[NF, Nc] = size(data);

[B, A] = butter(order, wn);


for I = 1:Nc,
    y = filtfilt(B, A, data(:,I));
    Y(:,I) = y;			 				
end
