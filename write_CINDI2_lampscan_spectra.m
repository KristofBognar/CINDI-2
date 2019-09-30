%% To write lamp scan spectra for Sebastian Donner's elevation calibration
% study

% select full scans manually, write one spectrum per file
% file structure similar to CINDI-2 raw spectra submission (Preprocs_CINDI_writespec)

%% Lamp scans have no altitude correction!! 
% cor_alt for normal measurements was -0.2
% for lamp scans, alt = alt + cor_alt = alt + 0.2 

%% Load calibration file (from QDOAS)
calibration=load('/home/kristof/work/GBS/PEARL-GBS/CINDI-2_2016/data/CINDI2_spectra/cal_P1_2016_253_2.asc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% General input -- need to be in folder for given day

yyyy=2016;
mm=09;
dd=10; % day is either the 10th or the 19th

% make list of all files
temp = dir('*.mat'); 
f_list = {temp.name}; % cell array of file names

%% Read in all data
all_spec=[];
for i=1:length(f_list)
    
    load(f_list{i});
    all_spec=[all_spec; spec];
    
end

%% Write one ascii file per spectrum
for i=1:size(all_spec,1)
    
    %% Open output file
    filename=['UToronto_36_VIS_' num2str(yyyy) num2str(mm) num2str(dd) ...
              '_' num2str(i) '_v1.asc'];

    fid = fopen(filename, 'w');

    %% print file header
    fprintf(fid, '%s\n', '# Station = Cabauw (51.97 N,4.93 E)');
    fprintf(fid, '%s\n', '# Institute = University of Toronto');
    fprintf(fid, '%s\n', '# PI name = Kimberly Strong (strong@atmosp.physics.utoronto.ca)');
    fprintf(fid, '%s\n', '# Instrument = VIS');
    fprintf(fid, '%s\n', '# Size of the detector = 2048');

    %% print spectra

    % spectrum to be written
    tmp=all_spec(i,:);
    
    % start and end time
    start_time=tmp(1:3);
    end_time=tmp(4:6);

    % mean time
    mt = ((tmp(1)/24 + tmp(2)/(60*24) + tmp(3)/(3600*24)) +...
                     (tmp(4)/24 + tmp(5)/(60*24) + tmp(6)/(3600*24)))/2;
    mean_time=zeros(1,3);
    mean_time(1)=floor(mt*24);
    mean_time(2)=floor((mt-mean_time(1)/24)*60*24);
    mean_time(3)=(mt-mean_time(1)/24-mean_time(2)/(60*24))*(3600*24);
    
    % viewing elevation and azimuth
    v_el=tmp(7)+0.2;
    v_az=tmp(8);

    % total time (in seconds)
    start_str=[num2str(start_time(1)) ':' num2str(start_time(2)) ':' ...
               num2str(round(start_time(3)))];
    end_str=[num2str(end_time(1)) ':' num2str(end_time(2)) ':' ...
               num2str(round(end_time(3)))];
    total_time=(datenum(end_str,'HH:MM:SS')-datenum(start_str,'HH:MM:SS'))*24*3600;

    % integration time (in seconds) and n.o. scans
    no_ac=tmp(10);
    t1int=tmp(9)/1000;
    total_exp=t1int*no_ac;

    %% write spectrum header
    fprintf(fid,'%s\n', ['Date (DD/MM/YYYY) = ' ...
        sprintf('%02d/%02d/%04d', dd, mm, yyyy)]);

    fprintf(fid,'%s\n', ['UTC Time (hh:mm:ss) = ' ...
        sprintf('%02d:%02d:%02d', mean_time(1),mean_time(2),round(mean_time(3)))]);

    fprintf(fid,'%s\n', ['UTC Start Time (hh:mm:ss) = ' ...
        sprintf('%02d:%02d:%02d', start_time(1),start_time(2),round(start_time(3)))]);

    fprintf(fid,'%s\n', ['UTC End Time (hh:mm:ss) = ' ...
        sprintf('%02d:%02d:%02d', end_time(1),end_time(2),round(end_time(3)))]);

    fprintf(fid, '%s\n', ['Viewing Elevation Angle (deg) = ' num2str(v_el)]);
    fprintf(fid, '%s\n', ['Viewing Azimuth Angle (deg) = ' num2str(v_az)]);
    fprintf(fid, '%s\n', ['Total Measurement Time (sec) = ' num2str(total_time)]);
    fprintf(fid, '%s\n', ['Total Acquisition Time (sec) = ' num2str(total_exp)]);
    fprintf(fid, '%s\n', ['Exposure time (sec) = ' num2str(t1int)]);

    %% write wavelength calibration and spectrum
    out=[calibration'; tmp(11:end)];
    fprintf(fid, '%.6f \t %.6f \n', out);


    fclose(fid);
    
end