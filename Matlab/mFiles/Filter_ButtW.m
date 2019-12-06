function [markerStruct, forceStruct] = Filter_ButtW(markerStruct,mFc,mFs, forceStruct,fFc,fFs)
%Example: [markerStruct, forceStruct] = Filter_ButtW(markerStruct,mFc=15,mFs=100, forceStruct,fFc=15,fFs=2000)
%filter markerData
fc = mFc;
fs = mFs;
[b,a] = butter(4,fc/(fs/2));
fields = fieldnames(markerStruct);
for i=1:numel(fields)-1
%     Q(isnan(Q))=0;
buff = markerStruct.(fields{i});
buff(isnan(buff))=0;
  markerStruct.(fields{i}) = filter(b,a,buff);
end

%filter forceData
fc = fFc;
fs = fFs;
[b,a] = butter(4,fc/(fs/2));
fields = fieldnames(forceStruct);
for i=1:numel(fields)-1
    buff = forceStruct.(fields{i});
buff(isnan(buff))=0;

  forceStruct.(fields{i}) = filter(b,a,buff);
end
