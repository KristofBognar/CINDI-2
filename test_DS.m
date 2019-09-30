%% Read in the QDOAS data

cd /home/kristof/work/GBS/PEARL-GBS/CINDI-2_2016/QDOAS

load('visible_cols.mat')

% read in QDOAS output file -- have to name that file manually!
file_nm= ['ds_.ASC'];
fid = fopen(file_nm, 'r');
fgetl(fid);
qdoas_raw = (fscanf(fid,'%f', [tot_nbr,inf]))';
fclose(fid);

% Now sort the data up by fractional day and ID whether there are
% doubles of some values
qdoas_raw = sortrows(qdoas_raw, no2_vis.fd);
all_ind = 1:length(qdoas_raw(:,1));
[a, unique_ind,b] = unique(qdoas_raw(:, no2_vis.fd));
diff_ind = setdiff(all_ind, unique_ind);
if ~isempty(diff_ind)
    disp('[WARNING]: File contains multiple entries taken at the same time.')
    disp('frac_day')
    disp('---------')
    disp(qdoas_raw(diff_ind, no2_vis.fd))
end

% filter out rms = 0, and dscd = 9999
ind = find(...
    qdoas_raw(:, no2_vis.rms) ~= 0 & ...
    qdoas_raw(:, no2_vis.no2_dscd) ~= 9999 & ...
    qdoas_raw(:, no2_vis.no2_dscd_err) ~= 9999 & ...
    qdoas_raw(:, no2_alt.rms) ~= 0 & ...
    qdoas_raw(:, no2_alt.no2_dscd) ~= 9999 & ...
    qdoas_raw(:, no2_alt.no2_dscd_err) ~= 9999 & ...
    qdoas_raw(:, o3_vis.rms) ~= 0 & ...
    qdoas_raw(:, o3_vis.o3_dscd) ~= 9999 & ...
    qdoas_raw(:, o3_vis.o3_dscd_err) ~= 9999 );
qdoas_raw = qdoas_raw(ind,:);

