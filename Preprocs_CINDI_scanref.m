function Preprocs_CINDI_scanref(yyyy, jjj_beg, jjj_end, instrument, loc, ...
    version, measurement_mode, view_elev, view_azimuth, filter_id,... 
    sat_lev, nbr_pix, date_format)
%--------------------------------------------------------------------------
% Preprocs program modified for CINDI-2 requirements by Kristof Bognar
% Use this for writing .spe input for sequential ref spectrum analysis
%
%%% ONLY WRITES .spe FILES!!!
%
% 
% Preprocs_v6(yyyy, jjj_beg, jjj_end, instrument, loc, ...
%     version, measurement_mode, view_elev, view_azimuth, filter_id,... 
%     sat_lev, nbr_pix, date_format)
% 
% Preprocs.m
% Written by Joseph Mendonca April 2010
% Modified by Cristen Adams May 2010
% Modified by Xiaoyi Zhao June 2012, add new output file format, ASCII (.spe), for
% QDOAS. The .spe file could be use for auto-scan MAX-DOAS analysis.
% Modified by Kristof Bognar 
%   April 2016: extract ZS measurements from MAX-DOAS scans,
%       auto-replace some known bad el/az values, set el=90 for ZS mode, 
%   September 2016: add viewing azimuth to .spe files for 2D MAX-DOAS
%       measurements
%   September 2016: add new measurement modes for horizon scan and
%       almucantar scan. Only write one file for these (use write_files_hs_as function)
%
% This matlab script reads the data taken the julian days specified
% (according to UTC-time).  In then splits the data up into spectrum files
% according to solar midnight.
%
% The following directory structure/files are required, where * is your
% working directory:
% */csv : folder containing csv files
% */dark_c.inp : dark current file
% */bias.inp : bias file
%
% Note the preprocs.m file can be kept anywhere (ie: should not be copied
% to each individual work directory).  Just use the "addpath" command so
% that matlab can find the function.
%
% Output files will go to */Processed
%
% The required input variables are:
% - yyyy: integer year (ie: 2007)
% - jjj_beg: integer start julian day
% - jjj_end: integer end julian day
% - instrument: 'U'-is for UT-GBS , 'P'- is for PEARL-GBS
% - location:   'C'=Cabauw, 'E' =Eureka, 'T'=Toronto, 'V' = Vanscoy)
%
% The optional input variables are:
% - version: version of CSV file (integer).  
%                   0 = multiple columns for spectrum, each corresponding with
%                       area of the CCD
%                   1 = 2004 - 2009  (single column for spectra)
% -                 2 = 2010- (day 78 for PGBS) - extra headers
%                   99 = 1999 EUreka data
% - measurement_mode (integer): 0 = dark current, 1 = zenith, 2 = MAX-DOAS, 
%                     3 = direct-sun, 4 = lunar 
% - view_elev, view_azimuth: viewing elevation and viewing azimuth
% (tracker)
% - filter_id: filter identification # (describe later)
% - sat_lev: counts for level of satuaration (default = 60000)
% - nbr_pix: number of pixels (default = 2048)
% - date_format:
%
%
% Example:  For measurments taken by the UT-GBS in Eureka from days 100-110
% in 2007 type in the command line.  Be sure to run the software from the
% working directory!!!!
% Preprocs_v5(2008, 100, 110, 'U', 'E');
% Preprocs_v5(2011, 063, 079, 'U', 'E', 2);
%--------------------------------------------------------------------------

%% set input parameters

%Default values when user does not input these values.
if nargin<7
   disp('Default settings being used')
   measurement_mode = 1;
   view_elev = -99.99;
   view_azimuth= -99.99;
   filter_id=0;
   if nargin < 6
       version = 1;
   end
   sat_lev = 60000;
end
if nargin < 12,
    nbr_pix = 2048;
    date_format = 'dd/mm/yyyy';
end

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
     noon=17;
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
     noon=11; %noon hour in UTC
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

%% loop over desired julian days
% now loop over the desired Julian days
for jjj = jjj_beg:jjj_end

    % get the month, day corresponding with current julian
    [dd, mm] = Julian2Date(yyyy, jjj);

    % loop over hours in a day and retrieve csv information into those days
    % this will fill the rows of the line3, line4, line5, and sp matrices
    for hh = 0:23,

        file_nm = csv_file_nm(yyyy, mm, dd, hh);
        [line3_tmp, line4_tmp, line5_tmp, line6_tmp, sp_tmp, is_sat_tmp, version, nbr_cols] =...
            read_csv(file_nm, location, version, sat_lev, nbr_pix, date_format, nbr_cols);

        % add this to the running line3 matrix 
        line3 = [line3; line3_tmp];
        line4 = [line4; line4_tmp];
        line5 = [line5; line5_tmp];
        line6 = [line6; line6_tmp];
        try
            sp = [sp; sp_tmp];
        catch
            disp('Error reading spectra')
            return
        end
        is_sat = [is_sat; is_sat_tmp];
        
    end
    
    % if there is no data for a given day, then skip to the next day
    if isempty(line4),
        disp(['No csv data for day: ' num2str(jjj)])
        continue
    end
    
    % Now get relevant indices for midnight etc.  

    % find rows within the midnight hour range on the current julian day
    % according to UTC
    L = length(line4(:,1));
    indtot = find( (line3(:,4) > (midnight_hr-1)) & ...
        (line3(:,4) < (midnight_hr+1)) & (line3(:,1) == dd) & ...
        (line3(:,2) == mm) & (line3(:,3) == yyyy));
    
    % If there are no rows matching these criteria, find rows with time
    % before the midnight hour on the current julian day according to UTC.
    % If there are no rows with this requirement, then call the first 
    % index midnight. 
    if isempty(indtot),
        % indices before midnight hour on current day (or from previous
        % day)
        ind_bef = find( (line3(:,4) < midnight_hr) | ...
        ( (datenum(fliplr(line3(:,1:3))) - datenum(yyyy, mm, dd)) == -1 ) );
        % If there is nothing before midnight time, skip to the next julian
        % day, keeping whatever data we have to be recorded on the next
        % day.
        if isempty(ind_bef)
            disp(['No data before midnight for day: ' num2str(jjj)])
            continue
        % Otherwise, the last of the values before midnight time-range
        % because the midnight index.
        else
            midnight_ind = ind_bef(length(ind_bef));
        end
        
    % If there are within the midnight range, then find the index of the
    % minimum in this range.
    else
        midnight_elev = min(line4(indtot, 1));
        midnight_ind = find(line4(:,1) == midnight_elev);
    end
    
    % these are our vectors for one day according to solar midnight
    line3_fin = line3(1:midnight_ind,:);
    line4_fin = line4(1:midnight_ind,:);
    line5_fin = line5(1:midnight_ind,:);
    if version == 2
        line6_fin = line6(1:midnight_ind,:);
    else
        line6_fin=[line5_fin(:,1),...
            ones(length(line5_fin(:,1)),1) * [view_elev, view_azimuth, filter_id]]; % make line 6
    end
    line7=[location.longitude, location.latitude];
    is_sat_fin = is_sat(1:midnight_ind,:);
    sp_l1 = sp(1:midnight_ind,:);
    
    % calculate the l2 spectrum by subtracting the DC and the bias
    sp_l2 = calc_l2(sp(1:midnight_ind,:), line5(8), nbr_pix);
    
    % If midnight hour < 12UTC, name the file by the previous day
    if midnight_hr < 12, 
        jjj_str = num2str(jjj-1);
    else
        jjj_str = num2str(jjj);
    end
    if length(jjj_str) == 2, jjj_str = ['0' jjj_str]; end
    if length(jjj_str) == 1, jjj_str = ['00' jjj_str]; end
    
    % Make the dark current files by finding the dark current indices and
    % then using the write_files function
    ind = find(line5_fin(:,1) == 0); % closed shutter
%     if ~isempty(ind),
%         filename = ...
%             [instrument loc '0D_' num2str(yyyy) '_' jjj_str];
%         write_files(filename, ind, line3_fin, line4_fin, line5_fin,...
%             line6_fin, line7, sp_l1, sp_l2);
%     end

    % line5 columns: shutter, ideal no counts, slit, groove, turret, blaze,
    % centre, integration time, no accum, mean, min, max TCCD, TBOX
    
    % Kristof: Replace some known bogus elevation / azimuth values
    % line6 columns: meas. mode, viewing elev, viewing az, filter
    bad_ind_el = find(line6_fin(:,2) == 9000000.00 | line6_fin(:,2) == 90000000.00);
    if ~isempty(bad_ind_el), line6_fin(bad_ind_el,2) = 90.00; end

    bad_ind_az = find(line6_fin(:,3) == 3500000.00 | line6_fin(:,3) == 35000000.00);
    if ~isempty(bad_ind_az), line6_fin(bad_ind_az,3) = 35.00; end

    
    % Make the solar files by finding the open shutter indices and correct
    % turret indices, then use the write_files function
    for turr = 0:2,
%         mode_vec = ['Z'; 'O'; 'S'; 'M'; 'X'];
%         for mode = 1:5
        mode_vec = ['Z'; 'O'; 'S'; 'M'; 'X'; 'H'; 'A'];
        for mode = 1:7
            if mode ~= 5
                ind = find(line6_fin(:,1) == mode & line5_fin(:,5) == turr & ...
                    line5_fin(:,1) == 1 & is_sat_fin == false);

                % include almucantar scans and horizon scans in MAXDOAS file                
                if mode == 2 
                    indhs = find(line6_fin(:,1) == 6 & line5_fin(:,5) == turr & ...
                        line5_fin(:,1) == 1 & is_sat_fin == false);
                    if ~isempty(indhs), ind=sort([ind;indhs]); end
                    indas = find(line6_fin(:,1) == 7 & line5_fin(:,5) == turr & ...
                        line5_fin(:,1) == 1 & is_sat_fin == false);
                    if ~isempty(indas), ind=sort([ind;indas]); end
                end
   
                if mode == 1
                    % replace all elevation values for ZS
                    % measurements in case labview recorded strange numbers
                    line6_fin(ind,2) = 90.00;
                    
                    % find zenith sky measurements in max doas
                    % scans and add them to Z files
                    ind_zenith = find(line6_fin(:,1) == 2 & ...
                        line5_fin(:,5) == turr & line5_fin(:,1) == 1 & ...
                        is_sat_fin == false & line6_fin(:,2) == 90.00);
                    
                    if ~isempty(ind_zenith), ind = sort([ind;ind_zenith]); end
                    
                    if loc=='E' && ~isempty(ind)
                        % average noon spectra
                       [ind,  line3_fin, line4_fin, line5_fin, line6_fin, ...
                            sp_l1, sp_l2] = avg_noon(ind, line3_fin, line4_fin, ...
                            line5_fin, line6_fin, sp_l1, sp_l2, noon,turr,is_sat_fin,loc);
                    
                        % average twilight spectra
                       [ind,  line3_fin, line4_fin, line5_fin, line6_fin, ...
                            sp_l1, sp_l2] = avg_twilight(ind, line3_fin, line4_fin, ...
                            line5_fin, line6_fin, sp_l1, sp_l2, noon,turr,is_sat_fin,loc);
                    
                    end
                end
                
                % set to true if averaging of noon reference is required
                % mode and time range are hard coded
%                 if mode == 2 && loc=='C' && ~isempty(ind)
%                     % put all zs measurements in the maxdoas file
%                     ind_zenith = find(line6_fin(:,1) == 1 & line5_fin(:,5) == turr & ...
%                         line5_fin(:,1) == 1 & is_sat_fin == false);
%                     if ~isempty(ind_zenith), ind=sort([ind;ind_zenith]); end
%                     
%                     % average noon spectra and reset index accordingly
%                     [ind,  line3_fin, line4_fin, line5_fin, line6_fin, ...
%                     sp_l1, sp_l2] = avg_noon(ind, line3_fin, line4_fin, ...
%                     line5_fin, line6_fin, sp_l1, sp_l2, noon,turr,is_sat_fin,loc);
%                 end

                
            elseif mode == 5
                ind = find(line5_fin(:,5) == turr & is_sat_fin == true);
            end
            
%             if mode < 6
            if mode == 2 % only write maxdoas files
                if ~isempty(ind),
                    filename = ...
                        [instrument loc num2str(turr) mode_vec(mode) ...
                        '_' num2str(yyyy) '_' jjj_str];
                    write_files(filename, ind, line3_fin, line4_fin, line5_fin,...
                        line6_fin, line7, sp_l1, sp_l2);
                end
%             elseif mode == 6
%                 if ~isempty(ind),
%                     filename = ...
%                         [instrument loc num2str(turr) mode_vec(mode) ...
%                         '_' num2str(yyyy) '_' jjj_str];
%                     write_files_hs_as(filename, ind, line3_fin, line4_fin, line5_fin,...
%                         line6_fin, line7, sp_l1, sp_l2);
%                 end
                
            end
        end
    end
    
    if length(sp(:,1)) ~= length(line3(:,1))
        disp('Error: something has gone terribly terribly wrong on this day!')
        return
    end
    
    
    % carry over any data that has not yet been written a file to the next
    % day
    if midnight_ind < L,
        line3 = line3(midnight_ind+1:L,:);
        line4 = line4(midnight_ind+1:L,:);
        line5 = line5(midnight_ind+1:L,:);
        if ~isempty(line6)
            line6 = line6(midnight_ind+1:L,:);
        end
        sp = sp(midnight_ind+1:L,:);
        is_sat = is_sat(midnight_ind+1:L,:);
    else
        line3 = [];
        line4 = [];
        line5 = [];
        line6 = [];
        sp = [];
    end
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
%     fid(1) = fopen([filename '.L1'], 'w');
%     fid(2) = fopen([filename '.L2'], 'w');
    fid(1) = fopen([filename '.spe'], 'w');%.spe file for auto-scan MAX-DOAS
    fid(2) = fopen([filename '.tint'], 'w');% file to save MAX-DOAS integration times and n.o. scans
    
    j = 0;
    for j = 1:length(ind)
        
        i = ind(j); % index of line matrices
        
        line3 = line3_mat(i,:);
        line4 = line4_mat(i,:);
        line5 = line5_mat(i,:);
        line6 = line6_mat(i,:);

        sza = 90 - line4(1,1); % calculate SZA for .spe file
        fraction_time = line3(1,4)+line3(1,5)/60+line3(1,6)/3600; %calculate fraction time for .spe file
        
        % old version: no viewing azimuth written in file
        % write frist 4 column [SZA, viewing_elevation, dd/mm/yyyy, fraction_time]
        % fprintf(fid(3), '%d %d %d/%d/%d %d ',sza,line6(1,2),line3(1,1),line3(1,2),line3(1,3),fraction_time);

        % Kristof: write viewing azimuth to file for 2D MAX-DOAS
        % measurements (CINDI-2)
        % write frist 5 column [SZA, viewing azimuth, viewing_elevation, dd/mm/yyyy, fraction_time]
        fprintf(fid(1), '%d %d %d %d/%d/%d %d ',sza,line6(1,3),line6(1,2),line3(1,1),line3(1,2),line3(1,3),fraction_time);
        % write tint and scans to file with frac time
        fprintf(fid(2), '%d %d %d %d ',j,fraction_time,line5(1,8),line5(1,9));
        fprintf(fid(2), '%d\n', []);
        
        % write spec to .spe
        fprintf(fid(1), '%.1f ', sp_l2(i,:));
        fprintf(fid(1), '%d\n', []);
    end

    fclose(fid(1));fclose(fid(2));

end
%%
function write_files_hs_as(filename, ind,  line3_mat, line4_mat, line5_mat, line6_mat,...
    line7, sp_l1, sp_l2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% same as write_files, but for horizon and almucantar scans, only saves
% .spe-like file, no L1 or L2
    
    disp(['Writing to ' filename])
    fid = [];
    fid(1) = fopen([filename '.hs'], 'w');
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
        
        % write spec to file
        fprintf(fid(1), '%.1f ', sp_l2(i,:));
        % write elev, az and time after spectrum
        fprintf(fid(1), '%d %d %d ',line6(1,2),line6(1,3),fraction_time);
        fprintf(fid(1), '%d\n', []);
    end

    fclose(fid(1));%fclose(fid(2));

end
%%
function [ind_out,  line3_out, line4_out, line5_out, line6_out, ...
    sp_l1_out, sp_l2_out] = avg_noon(ind, line3_fin, line4_fin, ...
    line5_fin, line6_fin, sp_l1, sp_l2, noon,turr,is_sat_fin,loc)

    % to find noon zs spectra and average them 
    
    % pre-assign output in case index is empty
    ind_out=ind;
    line3_out=line3_fin;
    line4_out=line4_fin;
    line5_out=line5_fin;
    line6_out=line6_fin;
    sp_l1_out=sp_l1;
    sp_l2_out=sp_l2;
    
    ind_noon=[];
    
    %find indices of ZS spectra around local noon
    if loc=='C'
        % avg between 11:30 and 11:40 as per CINDI-2 semiblind protocol
        ind_noon=find(line3_fin(:,10)==noon & line3_fin(:,11)>=30 & ...
            line3_fin(:,11)<=40 & line6_fin(:,1)==1 & ...
            line5_fin(:,5) == turr & line5_fin(:,1) == 1 & ...
            is_sat_fin == false & line6_fin(:,2) == 90);
    elseif loc=='E'
        % noon is around 17:40-17:50 UTC -- too variable and there's not
        % enough ZS measurements due to MAXDOAS scans: use min SZA and a
        % given range
        % NOTE: 90 degree measurements are part of MAX-DOAS (mode=2)
        
        % this finds largest solar elev (smallest SZA) among all spectra
        % (including MAXDOAS)
        % mean solar elevation for spectrum is line4(:,1)
        [noon_elev,noon_elev_ind]=max(line4_fin(:,1));
        
        % check if highest solar elev is near actual noon (in case of
        % missing data for the day)
        if line3_fin(noon_elev_ind,10)==noon && line3_fin(noon_elev_ind,11)>=30
            ind_noon=find(line4_fin(:,1) > noon_elev-0.03 & ...
                line6_fin(:,2)==90 & ...
                line5_fin(:,5) == turr & line5_fin(:,1) == 1 & ...
                is_sat_fin == false)
        end
    end
    
    % average noon spectra
    if ~isempty(ind_noon)
        % tint and scans should be added, not averaged
        tint=sum(line5_fin(ind_noon,8));
        no_ac=sum(line5_fin(ind_noon,9));

        % first start and last end times should be kept
        start_time=line3_fin(ind_noon(1),10:12);
        end_time=line3_fin(ind_noon(end),16:18);
        date=line3_fin(ind_noon(end),1:3);

        % recalculate average time
        t1=start_time(1)+start_time(2)/60 + start_time(3)/3600;
        t2=end_time(1)+end_time(2)/60 + end_time(3)/3600;
        a1=floor(mean([t1,t2]));
        a2=(mean([t1,t2])-a1)*60;
        a3=(a2-floor(a2))*60;
        avg_time=[a1,floor(a2),round(a3)];
            
%         elseif loc=='E'
%             % use values of smallest SZA measurement
%             start_time=line3_fin(noon_elev_ind,10:12);
%             end_time=line3_fin(noon_elev_ind,16:18);
%             avg_time=line3_fin(noon_elev_ind,4:6);
%             date=line3_fin(noon_elev_ind,1:3);
%         end
        
        % average all other parameters
        line4_t=mean(line4_fin(ind_noon,:),1);
        line5_t=mean(line5_fin(ind_noon,:),1);
        line6_t=mean(line6_fin(ind_noon,:),1);

        % average spectra
        sp_l1t=mean(sp_l1(ind_noon,:),1);
        sp_l2t=mean(sp_l2(ind_noon,:),1);

        % assign new values (same values on each line, so the indices of
        % rows after noon are notmessed up)
        for ijk=1:size(ind_noon,1)
            line3_fin(ind_noon(ijk),:)=[date,avg_time,date,start_time,date,end_time];
            line4_fin(ind_noon(ijk),:)=line4_t;
            line5_fin(ind_noon(ijk),:)=line5_t;
            line6_fin(ind_noon(ijk),:)=line6_t;

            sp_l1(ind_noon(ijk),:)=sp_l1t;
            sp_l2(ind_noon(ijk),:)=sp_l2t;
        end
        % add correct tint and scans values
        line5_fin(ind_noon,8)=tint;
        line5_fin(ind_noon,9)=no_ac;


        % since now all lines around noon are the same, add
        % only one to the main ind variable
        ind=setdiff(ind,ind_noon); % remove averaged indices first
        ind=sort([ind;ind_noon(1)]); % add one line back in
        
        % assign final results
        ind_out=ind;
        line3_out=line3_fin;
        line4_out=line4_fin;
        line5_out=line5_fin;
        line6_out=line6_fin;
        sp_l1_out=sp_l1;
        sp_l2_out=sp_l2;
        
    else disp('   no noon measurements -- no averaging')    
    end
end
%%
function [ind_out,  line3_out, line4_out, line5_out, line6_out, ...
    sp_l1_out, sp_l2_out] = avg_twilight(ind, line3_fin, line4_fin, ...
    line5_fin, line6_fin, sp_l1, sp_l2, noon,turr,is_sat_fin,loc)

    % to find twilight zs spectra and average them 
    
    % pre-assign output in case index is empty
    ind_out=ind;
    line3_out=line3_fin;
    line4_out=line4_fin;
    line5_out=line5_fin;
    line6_out=line6_fin;
    sp_l1_out=sp_l1;
    sp_l2_out=sp_l2;
    
    ind_twilight=[];
    
    % find noon index (or index of largest solar elevation -- smallest SZA)
    [~,noon_elev_ind]=max(line4_fin(:,1));
    
    binsize=0.5; % width of SZA bins

    for ampm=0:1
    
        % separate morning and evening spectra by noon index
        % mean solar elevation for spectrum is line4(:,1)
        if ampm==0
            ind1=1;
            ind2=noon_elev_ind;
            min_elev=min(line4_fin(ind1:ind2,1));
        elseif ampm==1
            ind1=noon_elev_ind;
            ind2=ind(end);
            min_elev=min(line4_fin(ind1:ind2,1));
        end
                
        % loop over each bin and average spectra for SZA>=84
        if min_elev<6
            for bin_no=1:floor((6-min_elev)/binsize)+1
                
                % bin limits
                top=min_elev + bin_no*binsize;
                bottom=top-binsize;

                % find spectra in each bin
                ind_twilight=find(line4_fin(ind1:ind2,1)< top & line4_fin(ind1:ind2,1)>= bottom & ...
                    line6_fin(ind1:ind2,1)==1 & ...
                    line5_fin(ind1:ind2,5) == turr & line5_fin(ind1:ind2,1) == 1 & ...
                    is_sat_fin(ind1:ind2) == false);
                
                % make sure index refers to full array
                ind_twilight=ind_twilight+ind1-1;

                % average twilight spectra in each bin
                if ~isempty(ind_twilight)
                    % tint and scans should be added, not averaged
                    tint=sum(line5_fin(ind_twilight,8));
                    no_ac=sum(line5_fin(ind_twilight,9));

                    % first start and last end times should be kept
                    start_time=line3_fin(ind_twilight(1),10:12);
                    end_time=line3_fin(ind_twilight(end),16:18);
                    date=line3_fin(ind_twilight(end),1:3);

                    % recalculate average time
                    t1=start_time(1)+start_time(2)/60 + start_time(3)/3600;
                    t2=end_time(1)+end_time(2)/60 + end_time(3)/3600;
                    a1=floor(mean([t1,t2]));
                    a2=(mean([t1,t2])-a1)*60;
                    a3=(a2-floor(a2))*60;
                    avg_time=[a1,floor(a2),round(a3)];

                    % average all other parameters
                    line4_t=mean(line4_fin(ind_twilight,:),1);
                    line5_t=mean(line5_fin(ind_twilight,:),1);
                    line6_t=mean(line6_fin(ind_twilight,:),1);

                    % average spectra
                    sp_l1t=mean(sp_l1(ind_twilight,:),1);
                    sp_l2t=mean(sp_l2(ind_twilight,:),1);

                    % assign new values (same values on each line, so the indices of
                    % rows after noon are notmessed up)
                    for ijk=1:size(ind_twilight,1)
                        line3_fin(ind_twilight(ijk),:)=[date,avg_time,date,start_time,date,end_time];
                        line4_fin(ind_twilight(ijk),:)=line4_t;
                        line5_fin(ind_twilight(ijk),:)=line5_t;
                        line6_fin(ind_twilight(ijk),:)=line6_t;

                        sp_l1(ind_twilight(ijk),:)=sp_l1t;
                        sp_l2(ind_twilight(ijk),:)=sp_l2t;
                    end
                    % add correct tint and scans values
                    line5_fin(ind_twilight,8)=tint;
                    line5_fin(ind_twilight,9)=no_ac;


                    % since now all lines are the same, add
                    % only one to the main ind variable
                    ind=setdiff(ind,ind_twilight); % remove averaged indices first
                    ind=sort([ind;ind_twilight(1)]); % add one line back in

                    % assign final results
                    ind_out=ind;
                    line3_out=line3_fin;
                    line4_out=line4_fin;
                    line5_out=line5_fin;
                    line6_out=line6_fin;
                    sp_l1_out=sp_l1;
                    sp_l2_out=sp_l2;
                end    

            end
        end
    end
end
