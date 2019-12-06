%% Calculte MoS 
% I think I need to use trc and mot files from after exporting to OpenSim
% to be consistante in the idea of c3d to OpenSim...

fields = fieldnames(markerStruct);
time = markerStruct.('time');
% Get marker data using pelvic markers
for i = 1:4 % I think this is not right it start from 4:7
buf(:,:,i) = markerStruct.(fields{i});
end
com = mean(buf,3);

%% Calculate pendulum's length and determine eigenfrequency
% l = mean(tpose(28:31, 2));
% max Y direction of com roughly I need to correct this
l = max(com(:,2))/1000;
inv_pend_eigenfreq = sqrt(l / 9.81);

%% Calculate XCoM (for x and z directions)
com_vel_x = diff(com(:, 1)/1000)./diff(markerStruct.(fields{end}));
xcom_x = com(1:end-1, 1) + (com_vel_x * inv_pend_eigenfreq);

com_z_tm = com(:, 3) * -1; 
com_vel_z_tm = diff(com(:, 3)/1000)./diff(markerStruct.(fields{end}));%differentiate_data(com(:, 3));
xcom_z_tm = com_z_tm(1:end-1, 1)/1000 + (com_vel_z_tm * inv_pend_eigenfreq);

belt_dist = 1.5; %get this from DFLOW later, I just assumeed a speed for simplicity
com(:, 3) = com(:, 3) * -1 + belt_dist; % convert to "overground" values
com_vel_z = diff(com(:, 3)/1000)./diff(markerStruct.(fields{end})); %differentiate_data(com(:, 3));
xcom_z = com(1:end-1, 3)/1000 + (com_vel_z * inv_pend_eigenfreq);

%% z: forward, x: to left, y: up
mos_AP = zeros(1,length(markerStruct.(fields{end})(1:end-1)));
mos_ML = zeros(1,length(markerStruct.(fields{end})(1:end-1)));
for i = 1:length(markerStruct.('LMT2'))-1
if markerStruct.('LMT2')(i,3) > markerStruct.('RMT2')(i,3) % MoS
    mos_AP(1,i) =  markerStruct.('LMT2')(i,3)/1000 - xcom_z(i);
    mos_ML(1,i) =  markerStruct.('LLM')(i,3)/1000 - xcom_z(i);
else
    mos_AP(1,i) =  markerStruct.('RMT2')(i,3)/1000 - xcom_z(i);   
    mos_ML(1,i) =  markerStruct.('RLM')(i,3)/1000 - xcom_z(i);   
end 
end

%% plotting mos
limit = [0 1];
subplot(1,2,1); plot(time(1:end-1),mos_AP);title('AP MoS');ylim(limit);xlabel('Time (s)')
subplot(1,2,2); plot(time(1:end-1),mos_ML);title('ML MoS');ylim(limit);xlabel('Time (s)')

figure % plot LMT2 
plot(time,markerStruct.('LMT2')(:,3)/1000) 
