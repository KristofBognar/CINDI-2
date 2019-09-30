function average_noon(yyyy, jjj, instrument, loc, version)
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

%% set input parameters

filter_id=0;
sat_lev = 60000;

nbr_pix = 2048;
date_format = 'dd/mm/yyyy';

% Set the location parameters
% midnight_hr is the appromate time of solar midnight in UTC (within +/- 1
% hour
if loc=='T'
     location.longitude = -79.40;
     location.latitude = 43.66;
     location.altitude = 174;
     midnight_hr = 5;
elseif loc=='E'
     location.longitude = -86.416;
     location.latitude = 80.053;
     location.altitude = 610;
     midnight_hr = 5;
elseif loc=='V'
     location.longitude = -106.98;
     location.latitude = 52;
     location.altitude = 516;
     midnight_hr = 4; % ??
elseif loc=='C'
     location.longitude = 4.898;
     location.latitude = 51.96;
     location.altitude = 200;
     midnight_hr = 0;
elseif loc=='R'
    location.longitude = -94.87;
    location.latitude = 74.68;
    location.altitude = 0;
    midnight_hr = 5;
else
    disp('Error: location not recognized')
end     

% make output file directory if necessary
% data_dir = 'Processed/';
% [s,mess,messid] = mkdir('.',data_dir);

% Matrices that will be filled with L1 and L2 fields until files are written
line3 = [];
line4 = [];
line5 = [];
line6 = [];
sp = [];
is_sat = [];

    
nbr_cols = 4;

%%
[dd, mm] = Julian2Date(yyyy, jjj);

if loc=='C'
    hh = 11;
else
    error('Please set noon time');
end
    
file_nm = csv_file_nm(yyyy, mm, dd, hh);
[line3, line4, line5, line6, sp_tmp, is_sat, version, nbr_cols] =...
    read_csv(file_nm, location, version, sat_lev, nbr_pix, date_format, nbr_cols);

try
    sp = [sp; sp_tmp];
catch
    disp('Error reading spectra')
    return
end

        
   
    
    % these are our vectors for one day according to solar midnight
    line3_fin = line3;
    line4_fin = line4;
    line5_fin = line5;
    line6_fin = line6;
    line7=[location.longitude, location.latitude];
    is_sat_fin = is_sat;
    sp_l1 = sp;
    
    % calculate the l2 spectrum by subtracting the DC and the bias
    sp_l2 = calc_l2(sp, line5(8), nbr_pix);
             
    % line5 columns: shutter, ideal no counts, slit, groove, turret, blaze,
    % centre, integration time, no accum, mean, min, max TCCD, TBOX
    
    % Kristof: Replace some known bogus elevation / azimuth values
    % line6 columns: meas. mode, viewing elev, viewing az, filter
    bad_ind_el = find(line6_fin(:,2) == 9000000.00 | line6_fin(:,2) == 90000000.00);
    if ~isempty(bad_ind_el), line6_fin(bad_ind_el,2) = 90.00; end

    
    % average noon spectra
%     for turr=0:2
    for turr=1
        ind=find(line3_fin(:,11)>=30 & line3_fin(:,11)<=40 & line6_fin(:,1)==1 & ...
            line5_fin(:,5) == turr & line5_fin(:,1) == 1 & is_sat_fin == false & ...
            line6_fin(:,2) == 90);

        tint=sum(line5_fin(ind,8));
        no_ac=sum(line5_fin(ind,9));
        
        line3_fin=mean(line3_fin(ind,:),1);
        line4_fin=mean(line4_fin(ind,:),1);
        line5_fin=mean(line5_fin(ind,:),1);
        line6_fin=mean(line6_fin(ind,:),1);
        line7_fin=line7;
        
        line5_fin(8)=tint;
        line5_fin(9)=no_ac;
        
        sp_l1=mean(sp_l1(ind,:),1);
        sp_l2=mean(sp_l2(ind,:),1);
        
        if ~isempty(ind),
        filename = ...
            [instrument loc num2str(turr) '_noon_' num2str(yyyy) '_' num2str(jjj)];
        write_files(filename, 1, line3_fin, line4_fin, line5_fin,...
            line6_fin, line7, sp_l1, sp_l2);
        end
    end
   
      
    if length(sp(:,1)) ~= length(line3(:,1))
        disp('Error: something has gone terribly terribly wrong on this day!')
        return
    end
    
end

%% sub-functions
function [line3, line4, line5, line6, sp, is_sat, version, nbr_cols] = ....
    read_csv(filename, location, version, sat_lev, nbr_pix, date_format, nbr_cols)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reads all the spectra in a given csv file into lines that
% will correspond with the lines to be printed in the L1 and L2 files.
% This calculates the elevation angles for each given spectrum.
%
% Input:  filename (string), location (structure to be read in by the get
% sza functions.
%
% Output: line3, line4, line4, sp (all matrices with rows corresponding
% with individual spectra and columsn corresponding with entries.

    % init the vectors
    line3 = [];
    line4 = [];
    line5 = [];
    line6 = [];
    sp = [];
    is_sat = [];

    %opens the file and assigning a value to fid and check to see if a file
    %is found.  If not, continue to the next hour.
    fid = fopen(filename,'r');
    if fid == -1 %check to see if the file is found
        %disp(['No file: ' filename])
        return
    end
    
    % display filename
    disp(['Reading: ' filename])
    
    %skips first two headers
    fgetl(fid); fgetl(fid);
    
    line_nbr = 2;

    
    nbr_catches = 0;

    % now read through the file until we reach the end of the file
    while ~feof(fid)

        % read in header, remove ", and split it up into substrings by comma
        header_str = fgetl(fid);
        header_str = regexprep(header_str, '"', '');
        header = regexp(header_str, ',', 'split');
        
        if isempty(cell2mat(header)), continue, end
        
        if feof(fid), continue, end
        
        % read in the spectrum
        if version > 0 && version < 90
            sp_tmp = fscanf(fid, '%f', nbr_pix)';
            % and read in last empty line
            fgetl(fid);
            indy = find(sp_tmp > sat_lev);
            if length(indy) > 3;
                is_sat = [is_sat; true];
            else
                is_sat = [is_sat; false];
            end
            delta_line_nbr = nbr_pix+1;
        elseif version == 99
            sp_tmp_str = [];
            for i = 1:10
                sp_tmp_str = [sp_tmp_str ',' fgetl(fid)];
            end
            cell_tmp = regexp(sp_tmp_str,',','split');
            for w = 2:length(cell_tmp)
                sp_tmp(w) = str2num(cell_tmp{w});
            end
            nbr_pix = length(sp_tmp);
            indy = find(sp_tmp > sat_lev);
            if length(indy) > 3;
                is_sat = [is_sat; true];
            else
                is_sat = [is_sat; false];
            end
            delta_line_nbr = 11;
        elseif version == 0
            file_format = '%f';
            for i = 2:nbr_cols
                file_format = [file_format ',%f'];
            end
            sp_tmp=fscanf(fid, file_format, [nbr_cols,nbr_pix]);    
            [indx, indy] = find(sp_tmp > sat_lev);
            indx = unique(indx);
            sp_tmp(indx,:) = 0;
            sp_tmp = sum(sp_tmp);
            if length(indy > 3)
                is_sat = [is_sat; true];
            else
                is_sat = [is_sat; false];
            end
            delta_line_nbr = nbr_pix+1;
        end
        
        if isempty(sp_tmp); continue; end
        
        format_flag = 0;
        if length(sp_tmp) == 1
            format_flag = 1;
        elseif length(sp_tmp) < nbr_pix
            format_flag = 2;
        else
            sp_tmp = [sp_tmp, zeros(1,2048 - nbr_pix)];
            try
                st_date = cell2mat(header(1));
                if length(header) < 17
                    end_date = st_date;
                else
                    end_date = cell2mat(header(3));
                end
                if length(st_date) < 10
                    if version == 99
                        st_date = [st_date(1:6) '19' st_date(7:8)];
                        end_date = [end_date(1:6) '19' end_date(7:8)];
                    else
                        st_date =[st_date(1:6) '20' st_date(7:8)];
                        end_date =[end_date(1:6) '20' end_date(7:8)];
                    end
                    st_date_vec = [str2num(st_date(7:10)) ...
                        str2num(st_date(4:5)) str2num(st_date(1:2)) 0 0 0]; 
                    end_date_vec = [str2num(end_date(7:10)) ...
                        str2num(end_date(4:5)) str2num(end_date(1:2)) 0 0 0]; 
                else
                    st_date_vec = datevec(st_date, date_format);
                    end_date_vec = datevec(end_date, date_format);
                end
          
                sp = [sp; sp_tmp];
            catch
                format_flag = 1;
            end
        end
        
        flag_v0 = 0;
        if format_flag > 0
            nbr_catches = nbr_catches + 1;
            
            if nbr_catches > 4,
                disp('PROBLEM WITH FILE!')
                fclose(fid);
                return
            end
            
            disp('Encountered different file format')
            
            if nbr_catches == 2
                disp('Resetting number of spectral columns')
                fgetl(fid);
                sp_line = fgetl(fid);
                ind = findstr(',',sp_line);
                nbr_cols = length(ind)+1;
                flag_v0 = 1;
            end
            
            fclose(fid);
            fid = fopen(filename,'r');
            for i = 1:line_nbr; fgetl(fid); end
            if version == 0,
                if flag_v0 == 0
                    version = 1;
                end
                continue
            elseif version == 1,
                version = 0;
                continue
            else
                break
            end
        end
        
        nbr_catches = 0;
       
        tmp = datevec(header(2));%separate time
        st_date_vec(4:6) = tmp(4:6);
        
        if length(header) < 17
            end_date_vec(4:6) = st_date_vec(4:6);
        else
            tmp = datevec(header(4));%separate time
            end_date_vec(4:6) = tmp(4:6);
        end

        % now calculate the start and end sza
        st_elev = get_sol_elev(st_date_vec, location);
        end_elev = get_sol_elev(end_date_vec, location);

        % now calculate the average start dates and start times
        [avg_date_vec, avg_elev] = calculate_averages(st_date_vec,...
            end_date_vec, st_elev, end_elev);

        % now make line vectors corresponding to future L1 header lines
        % note that the date is flipped compared to the standard date vec
        % ie: dd, mm, yyyy, hh, mn, ss
        line3 = [line3;...
            fliplr(avg_date_vec(1:3)), avg_date_vec(4:6)...
            fliplr(st_date_vec(1:3)), st_date_vec(4:6), ...
            fliplr(end_date_vec(1:3)), end_date_vec(4:6)];
        line4 = [line4; avg_elev, st_elev, end_elev];

        % read in the fifth line of the L1 header
        tmp = [];
        if version == 99 && length(header) < 17
            for i = 3:15, tmp = [tmp str2num(cell2mat(header(i)))];  end
        else
            for i = 5:17, tmp = [tmp str2num(cell2mat(header(i)))];  end
        end
        line5 = [line5; tmp]; % first index is shutter

        % make the 6th line
        if version == 2,
            tmp = [];
            for i = 18:23, tmp = [tmp str2num(cell2mat(header(i)))];  end
            line6 = [line6; tmp(6) tmp(1) tmp(2) tmp(3)]; % first index is meas. mode
        else
            line6 = [];
        end        
        line_nbr = line_nbr + delta_line_nbr;
    end
    fclose(fid);
end
%%
function [avg_d, avg_e] = calculate_averages(d1, d2, e1, e2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the average time and average solar elevation
% angle for
% INPUT:  d1, d2 (input date vectors [yyyy, mm, dd, hh, mm, ss])
%         e1, e2 (input elevation angles: float)
% OUTPUT: avg_d (average date vector)
%         avg_e (average elevation angle)
%         

    % fill in necessary variables for function
    n1 = datenum(d1);
    n2 = datenum(d2);    
    n_avg = (n1+n2)/2;
    
    avg_d = datevec(n_avg);
    avg_e = (e1+e2)/2;
end
%%
function elev = get_sol_elev(date_vec, location)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function sets up the time structure and calls the sun_position
% function to calculate the solar elevation angle
%
% INPUT: date_vec [yyyy, mm, yy, dd, mm, hh, ss]
%        location (structure)
%
% OUTPUT: elev (float elevation angle of sun)

time.year = date_vec(1);
time.month = date_vec(2);
time.day = date_vec(3);
time.hour =date_vec(4);
time.min = date_vec(5);
time.sec = date_vec(6);
time.UTC = 0;
sun = sun_position(time, location);
elev = 90 - sun.zenith;

end
%%
function file_nm = csv_file_nm(yyyy, mm, dd, hh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates a csv filename for the date
%
% INPUT: yyyy, mm, dd, hh (integers for year, month, day, and hour of
% measurements in UTC.
%
% OUTPUT: csv filename (string)
    % create year string (last two digits only)
    yy_str = num2str(yyyy);
    yy_str = yy_str(3:4);

    % create month strings
    if mm < 10, mm_str = ['0' num2str(mm)];
    else mm_str = num2str(mm); end

    % create day strings
    if dd < 10, dd_str = ['0' num2str(dd)];
    else dd_str = num2str(dd); end

    % create hour string
    if hh < 10, hh_str = ['0' num2str(hh)];
    else hh_str = num2str(hh); end

    % create csv filename
    file_nm = ['./csv/' dd_str mm_str yy_str hh_str '.csv'];

end
%%
function sp_l2 = calc_l2(sp_l1, exp, nbr_pix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the L2 spectrum by subtracting dark-current and
% bias from an input spectrum.
%
% INPUT: exp (exposure time, s) 
%        sp_l1 (matrix where rows correspond with individual spectra and
% columns correspond pixel #)
%
% OUTPUT: sp_l2 (same as sp_l1 except corrected)
%
% NOTE: this function also requires that there be dark_c.inp and bias.inp
% directories in the working directory.

try 
    dc = textread('dark_c.inp','%f');
    bias = textread('bias.inp','%f');
catch
    disp('Error: dark_c.inp and bias.inp needs to be in working directory')
    user_input = input('Press enter to continue');
end

sp_l2 = [];
for i = 1:length(sp_l1(:,1))
    sp_l2(i, 1:nbr_pix) = sp_l1(i,1:nbr_pix) - dc' * exp /1000 - bias';
    if nbr_pix < 2048, 
        sp_l2(i, (nbr_pix+1):2048) = zeros(1,2048-nbr_pix);
    end
end

end
%%
function write_files(filename, ind,  line3_mat, line4_mat, line5_mat, line6_mat,...
    line7, sp_l1, sp_l2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the L1 file for a given day.  
%
% INPUT: filename - string filename for given dataset
%
%        ind - row indices corresponding with desired spectra (generally
%        before this spectra are filtered for measurement type, etc)
%
%        line3_mat, line4_mat, line5_mat - each row corresponds with a
%        different spectrum
%
%        line6, line7 - vectors for line 6 and line 7
%
%        sp_l1 and sp_l2 - level 1 and level 2 spectra respectively with
%        each row corresponding with a different spectrum and columns
%        corresponding with pixels
%
% OUTPUT:  no output except written .L1, .L2 and .spe files        
%
    
    disp(['Writing to ' filename])
    fid = [];
    fid(1) = fopen([filename '.L1'], 'w');
    fid(2) = fopen([filename '.L2'], 'w');
    fid(3) = fopen([filename '.spe'], 'w');%.spe file for auto-scan MAX-DOAS
    fid(4) = fopen([filename '.tint'], 'w');% file to save MAX-DOAS integration times and n.o. scans
    
    
    % print out total number of spectra
    for k = 1:2, fprintf(fid(k), '%d\n', length(ind)); end

    j = 0;
    for j = 1:length(ind)
        
        i = ind(j); % index of line matrices
        
        for k = 1:2
            fprintf(fid(k),'%s\n', '*********************************************************');

            % write line 1: number of spectrum
            fprintf(fid(k), '%d\n', j);

            %write line 2: pixels
            fprintf(fid(k), '%d\t %d\n', [0, 2047]);

            %write line 3: average date/time
            line3 = line3_mat(i,:);
            fprintf(fid(k), '%2.0f\t %u\t %4.0d \t%2.0f \t%2.0f \t%2.0f \t%d \t%d \t%d \t%d \t%d \t%2.0f \t%d \t%d \t%d \t%d \t%d \t%2.0f\n', line3);

            %write line 4 to file
            line4 = line4_mat(i,:);
            fprintf(fid(k), '%3.1f\t %3.1f\t %3.1f\n', line4);

            %writes line 5 to file
            line5 = line5_mat(i,:);
            fprintf(fid(k), '%.0f\t %.0f\t %.1f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t% .0f\t %.0f\t %.0f\t %.0f\t %.0f\n', line5);

            %write line 6 to file
            line6 = line6_mat(i,:);
            fprintf(fid(k), '%2.1f\t%4.2f\t%4.2f\t%2.1f\n', line6);

            %writes line 7 to file
            fprintf(fid(k), '%3.1f\t%3.1f\n', line7);
        end
        sza = 90 - line4(1,1); % calculate SZA for .spe file
        fraction_time = line3(1,4)+line3(1,5)/60+line3(1,6)/3600; %calculate fraction time for .spe file
        
        % old version: no viewing azimuth written in file
        % write frist 4 column [SZA, viewing_elevation, dd/mm/yyyy, fraction_time]
        % fprintf(fid(3), '%d %d %d/%d/%d %d ',sza,line6(1,2),line3(1,1),line3(1,2),line3(1,3),fraction_time);

        % Kristof: write viewing azimuth to file for 2D MAX-DOAS
        % measurements (CINDI-2)
        % write frist 5 column [SZA, viewing azimuth, viewing_elevation, dd/mm/yyyy, fraction_time]
        fprintf(fid(3), '%d %d %d %d/%d/%d %d ',sza,line6(1,3),line6(1,2),line3(1,1),line3(1,2),line3(1,3),fraction_time);
        % write tint and scans to file with frac time
        fprintf(fid(4), '%d %d %d %d ',j,fraction_time,line5(1,8),line5(1,9));
        fprintf(fid(4), '%d\n', []);
        
        % write spec to file
        fprintf(fid(1), '%.1f\n', sp_l1(i,:));
        fprintf(fid(2), '%.1f\n', sp_l2(i,:));
        % write spec to .spe
        fprintf(fid(3), '%.1f ', sp_l2(i,:));
        fprintf(fid(3), '%d\n', []);
    end

    fclose(fid(1)); fclose(fid(2));fclose(fid(3));fclose(fid(4));

end
%%
function write_files_hs_as(filename, ind,  line3_mat, line4_mat, line5_mat, line6_mat,...
    line7, sp_l1, sp_l2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% same as write_files, but for horizon and almucantar scans, only saves
% .spe-like file, no L1 or L2
    
    disp(['Writing to ' filename])
    fid = [];
    fid(1) = fopen([filename '.spe'], 'w');%.spe file for auto-scan MAX-DOAS
%     fid(2) = fopen([filename '.tint'], 'w');% file to save MAX-DOAS integration times and n.o. scans
    
    
    j = 0;
    for j = 1:length(ind)
        
        i = ind(j); % index of line matrices
            line3 = line3_mat(i,:);
            line4 = line4_mat(i,:);
            line5 = line5_mat(i,:);
            line6 = line6_mat(i,:);
            
        
         sza = 90 - line4(1,1); % calculate SZA for .spe file
        fraction_time = line3(1,4)+line3(1,5)/60+line3(1,6)/3600; %calculate fraction time for .spe file
 
 
        % write tint and scans to file with frac time
%         fprintf(fid(2), '%d %d %d %d ',j,fraction_time,line5(1,8),line5(1,9));
%         fprintf(fid(2), '%d\n', []);
        
        % write spec to .spe
        fprintf(fid(1), '%.1f ', sp_l2(i,:));
        % write elev, az and time after spectrum
        fprintf(fid(1), '%d %d %d ',line6(1,2),line6(1,3),fraction_time);
        fprintf(fid(1), '%d\n', []);
    end

    fclose(fid(1));%fclose(fid(2));

end
