function CINDI_2_NO2_format_all()

version ='v3';
% species ='NO2_460nm';  
species ='NO2_360nm';  
% data_source='instrument36';
data_source='MedianDSCD';

% data_path_all=['/home/kristof/work/profile_retrievals/profile_results/CINDI-2/commonDSCD_VIS_',version,'/no2/'];
data_path_all=['/home/kristof/work/profile_retrievals/profile_results/CINDI-2/commonDSCD_UV_',version,'/no2/'];

cd(data_path_all);
tmp=dir();
list={tmp.name};
N = length(list) - 2;
for i = 1:1:N
    % for submission of our data -- not necessary for common dSCDs    
% %     if strcmp(list(i+2,:), '20160915') | strcmp(list(i+2,:), '20160918') | strcmp(list(i+2,:), '20160921') | strcmp(list(i+2,:), '20160924') 
% %         missing_1st_elev = 1;
% %     else
% %         missing_1st_elev = 0;
% %     end
    missing_1st_elev = 0;
    data_path = [data_path_all list{i+2} '/'];
    CINDI_2_tg_avk_output_format(1,version,species,data_path,data_source);
    CINDI_2_tg_dSCD_output_format(1,version, species,data_path,missing_1st_elev,data_source);
    CINDI_2_tg_pf_output_format(1, version, species, data_path,data_source);
end

