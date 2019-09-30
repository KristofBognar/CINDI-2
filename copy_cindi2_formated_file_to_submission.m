function copy_cindi2_formated_file_to_submission()
%% remove or transfer?
transfer=true;

% destination = '/home/kristof/work/profile_retrievals/profile_results/CINDI-2/commonDSCD_VIS_v3/submission/';
% % source = '/home/kristof/work/profile_retrievals/profile_results/CINDI-2/commonDSCD_VIS_v3/aerosol/';
% source = '/home/kristof/work/profile_retrievals/profile_results/CINDI-2/commonDSCD_VIS_v3/no2/';

destination = '/home/kristof/work/profile_retrievals/profile_results/CINDI-2/commonDSCD_UV_v3/submission/';
source = '/home/kristof/work/profile_retrievals/profile_results/CINDI-2/commonDSCD_UV_v3/aerosol/';
% source = '/home/kristof/work/profile_retrievals/profile_results/CINDI-2/commonDSCD_UV_v3/no2/';

cd(source);
tmp=dir();
list={tmp.name};
N = length(list);
for i =3:1:N
    if transfer
        source_location = [source list{i} '/cindi2/*.asc'];
        copyfile( source_location , destination);
    else
        rm([source list{i} '/cindi2/']);
    end
end

% % % if ~transfer, rm(destination); end

