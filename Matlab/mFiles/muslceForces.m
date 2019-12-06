function [muscleforce] = muslceForces(SOs, trial , muscles)

%% To get mucle forces from SO structure using muscle group e.g. hamstring,
% which trial and input the SOs then the muscle force will be output.
% By Hossein Mokhtarzadeh 
% Date: Oct 8 2019

%% Example

% muscles = {'rhams','lhams','lquad','rgluts','lgluts','rgast','lgast','rsole','lsole','rTAs','lTAs'};
% trial = 'S+T01';

% [muscleforce] = muslceForces(SOs, 'S+T01' , 'rhams');


%% fields = {'glut_med1_r';'glut_med2_r';'glut_med3_r';'glut_min1_r';'glut_min2_r';'glut_min3_r';...
%     'semimem_r';'semiten_r';'bifemlh_r';'bifemsh_r';'sar_r';'add_long_r';'add_brev_r';'add_mag1_r';...
%     'add_mag2_r';'add_mag3_r';'tfl_r';'pect_r';'grac_r';'glut_max1_r';'glut_max2_r';'glut_max3_r';...
%     'iliacus_r';'psoas_r';'quad_fem_r';'gem_r';'peri_r';'rect_fem_r';'vas_med_r';'vas_int_r';...
%     'vas_lat_r';'med_gas_r';'lat_gas_r';'soleus_r';'tib_post_r';'flex_dig_r';'flex_hal_r';...
%     'tib_ant_r';'per_brev_r';'per_long_r';'per_tert_r';'ext_dig_r';'ext_hal_r';'glut_med1_l';...
%     'glut_med2_l';'glut_med3_l';'glut_min1_l';'glut_min2_l';'glut_min3_l';'semimem_l';'semiten_l';...
%     'bifemlh_l';'bifemsh_l';'sar_l';'add_long_l';'add_brev_l';'add_mag1_l';'add_mag2_l';'add_mag3_l';...
%     'tfl_l';'pect_l';'grac_l';'glut_max1_l';'glut_max2_l';'glut_max3_l';'iliacus_l';'psoas_l';...
%     'quad_fem_l';'gem_l';'peri_l';'rect_fem_l';'vas_med_l';'vas_int_l';'vas_lat_l';'med_gas_l';...
%     'lat_gas_l';'soleus_l';'tib_post_l';'flex_dig_l';'flex_hal_l';'tib_ant_l';'per_brev_l';...
%     'per_long_l';'per_tert_l';'ext_dig_l';'ext_hal_l';'ercspn_r';'ercspn_l';'intobl_r';'intobl_l';...
%     'extobl_r';'extobl_l';'FX';'FY';'FZ';'MX';'MY';'MZ';'hip_flexion_r_reserve';'hip_adduction_r_reserve';...
%     'hip_rotation_r_reserve';'knee_angle_r_reserve';'ankle_angle_r_reserve';'subtalar_angle_r_reserve';...
%     'mtp_angle_r_reserve';'hip_flexion_l_reserve';'hip_adduction_l_reserve';'hip_rotation_l_reserve';...
%     'knee_angle_l_reserve';'ankle_angle_l_reserve';'subtalar_angle_l_reserve';'mtp_angle_l_reserve';...
%     'lumbar_extension_reserve';'lumbar_bending_reserve';'lumbar_rotation_reserve';...
%     'calcn_l_ExternalForce_1_Fx';'calcn_l_ExternalForce_1_Fy';'calcn_l_ExternalForce_1_Fz';...
%     'calcn_l_ExternalForce_1_px';'calcn_l_ExternalForce_1_py';'calcn_l_ExternalForce_1_pz';...
%     'calcn_l_ExternalForce_1_Tx';'calcn_l_ExternalForce_1_Ty';'calcn_l_ExternalForce_1_Tz';...
%     'calcn_r_ExternalForce_2_Fx';'calcn_r_ExternalForce_2_Fy';'calcn_r_ExternalForce_2_Fz';...
%     'calcn_r_ExternalForce_2_px';'calcn_r_ExternalForce_2_py';'calcn_r_ExternalForce_2_pz';...
%     'calcn_r_ExternalForce_2_Tx';'calcn_r_ExternalForce_2_Ty';'calcn_r_ExternalForce_2_Tz';'time'};

%% get the fields of SO
fields = fieldnames(SOs(6).data);
for i=1:length(SOs)
trials(i) = string(extractBefore(SOs(i).name,'_'));
end
% get Ham, Quad, Sol, Gas, Glut, TA, muslces, use indexing

% muscles = [rhams;lhams];

switch muscles
    case 'rhams'
        mm = {'semimem_r','semiten_r', 'bifemlh_r','bifemsh_r'};
    case 'lhams'
        mm = {'semimem_l','semiten_l', 'bifemlh_l','bifemsh_l'};
    case 'rquad'
        mm = {'rect_fem_r','vas_med_r','vas_int_r','vas_lat_r'};
    case 'lquad'
        mm = {'rect_fem_l','vas_med_l','vas_int_l','vas_lat_l'};
    case 'rgluts'
        mm = {'glut_med1_r', 'glut_med2_r', 'glut_med3_r', 'glut_min1_r',...
            'glut_min2_r', 'glut_min3_r'};
    case 'lgluts'
        mm = {'glut_med1_l', 'glut_med2_l', 'glut_med3_l', 'glut_min1_l',...
            'glut_min2_l', 'glut_min3_l'};
    case 'rgast'
        mm = {'med_gas_r', 'lat_gas_r'};
    case 'lgast'
        mm = {'med_gas_l', 'lat_gas_l'};
    case 'rsole'
        mm = {'soleus_r'};
    case 'lsole'
        mm = {'soleus_l'};
    case 'rTAs'
        mm = {'tib_ant_r'};
    case 'lTAs'
        mm = {'tib_ant_l'};
end

[indmm, ~]=trialIndexFinder(string(mm),fields);
try 
[ind, ~]=trialIndexFinder(string(trial),trials);
muscleforce = 0;

for i = 1:length(indmm)
    Buff(:,i) = SOs(ind).data.(string(mm(i)));
end
muscleforce = sum(Buff,2);
        
catch 
    'check this later'
    muscleforce = NaN;
end






