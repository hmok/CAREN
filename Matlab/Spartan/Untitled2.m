%% z: forward, x: to left, y: up
mos_AP = zeros(1,length(markerStruct.(fields{end})(1:end-1)));
for i = 1:length(markerStruct.('LMT2'))
if markerStruct.('LMT2')(i,3) > markerStruct.('RMT2')(i,3) % MoS
    mos_AP(1,i) =  markerStruct.('LMT2')(i,3)/1000 - xcom_z(i);
  else 
 mos_AP(1,i) =  markerStruct.('RMT2')(i,3)/1000 - xcom_z(i);
end 
end
plot(mos_AP)