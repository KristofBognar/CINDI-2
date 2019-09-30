function [  ] = write_template( f_in, f_out )
%generic_out sample output file for CINDI2
%   

    %% Write header 

    % open file
    fid = fopen([f_out], 'w');

    % write necessary info (specific to retrieved trace gas)
    fprintf(fid(1), '%s\n', '* NofHeaderlines: ??');
    fprintf(fid(1), '%s\n', '* NofColumns: ?? (if any info missing, put -999, even if it is the whole column) ');
    fprintf(fid(1), '%s\n', '* Instrument identifier: UToronto_MAXDOAS');
    fprintf(fid(1), '%s\n', '* Retrieval code: ??QDOAS (v2.110, June 2015)??');
    fprintf(fid(1), '%s\n', '* Created by: Kristof Bognar');
    fprintf(fid(1), '%s\n', '* Version: ??_v1');
    fprintf(fid(1), '%s\n', '* X-Axis (Col 1) = Day of year (DOY) 2016 (please start with 0.0 for January 1st, 0:00 UTC)');
    fprintf(fid(1), '%s\n', '* Y1-Axis (Col 2) = Time of day in hours (UTC)');
    fprintf(fid(1), '%s\n', '* Y2-Axis (Col 3) = Total Integration Time(s) ');
    fprintf(fid(1), '%s\n', '* Y3-Axis (Col 4) = Solar Zenith Angle (°)');
    fprintf(fid(1), '%s\n', '* Y4-Axis (Col 5) = Solar Azimuth Angle (°) North=0, East=90');
    fprintf(fid(1), '%s\n', '* Y5-Axis (Col 6) = Elevation Angle (°)');
    fprintf(fid(1), '%s\n', '* Y6-Axis (Col 7) = Viewing Angle (°) North=0, East=90');
    % ... all QDOAS products for given retrieval (species+error)

    fprintf(fid(1), '%s\n', '* Y??-Axis (Col ??) = Ring');
    fprintf(fid(1), '%s\n', '* Y??-Axis (Col ??) = Ring_Error');
    fprintf(fid(1), '%s\n', '* Y??-Axis (Col ??) = Fit RMS (in OD)');
    fprintf(fid(1), '%s\n', '* Y??-Axis (Col ??) = Spectrum shift (nm, against FRS reference)');
    fprintf(fid(1), '%s\n', '* Y??-Axis (Col ??) = Relative Intensity (counts/integration time @ ???nm)');
    fprintf(fid(1), '%s\n', '* Y??-Axis (Col ??) = Colour index: (??? nm / ??? nm)');
    fprintf(fid(1), '%s\n', '* Y??-Axis (Col ??) = Intensity offset with normalisation by I, I is the mean intensity in the spectral analysis windows, constant term');
    fprintf(fid(1), '%s\n', '* Y??-Axis (Col ??) = Intensity offset, linear term');
    % other details
    fprintf(fid(1), '%s\n', '* Fit settings: 1');
    fprintf(fid(1), '%s\n', '* Fitting Window: ???-??? nm');
    fprintf(fid(1), '%s\n', '* Polynomial: ?? (?? coefficients)');
    fprintf(fid(1), '%s\n', '* Offset: ?? order ');
    fprintf(fid(1), '%s\n', '* Calibration: Based on reference SAO solar spectra (Chance and Kurucz, 2010) --> sao2010_solref_air.dat');
    fprintf(fid(1), '%s\n', '* Wavelength adjustment: all spectra shifted and stretched against reference spectrum');
    fprintf(fid(1), '%s\n', '* Reference: sequential (zenith spectrum of the scan)');
    % ... all the cross-section files go here
   
    % ... write column headers here
    head=sprintf('*DOY\tUTC\tTint\tSZA\tSAA\tElev\tViewing_angle ...');
    fprintf(fid(1), '%s\n', head);
    
    %fprintf(fid(1), '%s\n', '');
    
    fclose(fid);

    %% Read in the QDOAS data
    fid = fopen(file_nm, 'r');
    fgetl(fid);
    qdoas_raw = (fscanf(fid,'%f', [col.tot_nbr,inf]))';
    fclose(fid);
    
    % Now sort the data up by fractional day and ID whether there are
    % doubles of some values
    qdoas_raw = sortrows(qdoas_raw, col.fd);
    all_ind = 1:length(qdoas_raw(:,1));
    [a, unique_ind,b] = unique(qdoas_raw(:, col.fd));
    diff_ind = setdiff(all_ind, unique_ind);
    if ~isempty(diff_ind)
        disp('[WARNING]: File contains multiple entries taken at the same time.')
        disp('frac_day')
        disp('---------')
        disp(qdoas_raw(diff_ind, col.fd))
    end
    
    % filter out rms = 0, and dscd = 9999
    ind = find(qdoas_raw(:, col.rms) ~= 0 & ...
        qdoas_raw(:, col.dscd) ~= 9999 & ...
        qdoas_raw(:, col.err) ~= 9999 & ...
        qdoas_raw(:, col.sza) < filt.sza_max);
    qdoas_raw = qdoas_raw(ind,:);
    
    %% Create data array from QDOAS output
    
    % create output array
    out_array=zeros(size(qdoas_raw,1),colnum);
    
    % calculate total integration time (tint)
    tint=qdoas_raw(:,col.scans).*qdoas_raw(:,col.t_1int);
    
    % calculate relative intensity (ri) at given lambda
    ri=qdoas_raw(:,col.I_)./qdoas_raw(col.t_1int);
    
    % calculate color index (ci, I_lowest/I_highest)
    ci=qdoas_raw(:,col.I_)./qdoas_raw(col.I_);
    
    % calculate intensity offset
    
    
    % Assign columns
    out_array(:,1)=qdoas_raw(:,col.day);
    out_array(:,1)=qdoas_raw(:,col.time);
    out_array(:,1)=tint;
    out_array(:,1)=qdoas_raw(:,col.sza);
    out_array(:,1)=qdoas_raw(:,col.saa);
    out_array(:,1)=qdoas_raw(:,col.elev);
    out_array(:,1)=qdoas_raw(:,col.azim);
    
    %out_array(:,1)=qdoas_raw(:,col.); all the dscds
    
    out_array(:,1)=qdoas_raw(:,col.rms);
    out_array(:,1)=qdoas_raw(:,col.shift);
    out_array(:,1)=ri;
    out_array(:,1)=ci;
%     out_array(:,1)=qdoas_raw(:,col.);
    
    %% Writre data in the file

    dlmwrite(f_out,out_array,'delimiter','\t','precision',4,'-append')
         
end

