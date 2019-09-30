%%% Program to write CINDI-2 data files for resubmission using scan
%%% references -- only works if date and time fields are not in the QDOAS
%%% file

cd /home/kristof/work/GBS/PEARL-GBS/CINDI-2_2016/QDOAS/

load('visible_cols_scanref.mat');

for asd=256:272
%% Read in QDOAS data

% define measureent day
% day=270;
day=asd;

% get corresponding date
[dd,mm] = Julian2Date(2016,day);
date=datetime(2016,mm,dd,'Format','yyyyMMdd');
date=datestr(date,'yyyymmdd');

disp(['***   Writing data for ' num2str(mm) '/' num2str(dd) '   ***'])

% read in QDOAS output file -- have to name that file manually!
file_nm= ['maxdoas_scanref_' num2str(day) '.ASC'];
fid = fopen(file_nm, 'r');
fgetl(fid);
qdoas_raw = (fscanf(fid,'%f', [tot_nbr,inf]))';
fclose(fid);

% check if actually reading in the right day
if day~=round(qdoas_raw(1,no2_vis.fd))
   error('QDOAS file does not match selected day!') 
end

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
    qdoas_raw(:, no2_alt.no2_dscd_err) ~= 9999);
qdoas_raw = qdoas_raw(ind,:);

%% Get t_int and intensities

% read tint and scans from my custom file (only for sequential
% references, since .spe files don't save that info)
tintfile=['/home/kristof/work/GBS/PEARL-GBS/CINDI-2_2016/data/scanref/PC1O_2016_' num2str(day) '.tint'];
tintarray=load(tintfile);

% use spectrum number to match to QDOAS data (QDOAS numbers the spectra)
inds=intersect(tintarray(:,1),qdoas_raw(:,no2_vis.spec_no));

% tint is ms in my files!!!
t_1int=tintarray(inds,3)./1000;
scans=tintarray(inds,4);
% calculate total integration time (tint) in seconds
tint=scans.*t_1int;

% calculate relative intensity (ri) at given lambda (use exposure for
% single measurement here!)
ri_440=qdoas_raw(:,I_440)./t_1int;
ri_500=qdoas_raw(:,I_500)./t_1int;

% calculate color index (ci, as defined in data protocol)
ci_1=qdoas_raw(:,I_425)./qdoas_raw(I_440);
ci_2=qdoas_raw(:,I_412)./qdoas_raw(I_440);
ci_3=qdoas_raw(:,I_440)./qdoas_raw(I_500);

%% Assign selection-specific values

for selection=[1,2]
    
    if selection==1 % NO2 and O4 vis
        
        % set array size and output file name
        colnum=24;
        f_out=['CINDI_files_scanref/UToronto_MAXDOAS_36_NO2vis_CINDI2_' date '_v3.asc'];

        % create output array
        out_array=zeros(size(qdoas_raw,1),colnum);

        % Assign columns: general info -- same for all species
%         out_array(:,1)=day;
        out_array(:,1)=qdoas_raw(:,no2_vis.fd)-1;
        out_array(:,2)=qdoas_raw(:,no2_vis.ft);
        out_array(:,3)=tint;
        out_array(:,4)=qdoas_raw(:,no2_vis.sza);
        out_array(:,5)=qdoas_raw(:,no2_vis.saa)+180.0;
        out_array(:,6)=qdoas_raw(:,no2_vis.elev);
        out_array(:,7)=qdoas_raw(:,no2_vis.az);

        % DSCDs and errors, divided by factors as specified in header text
        out_array(:,8)=qdoas_raw(:,no2_vis.no2_dscd)./1e15;
        out_array(:,9)=qdoas_raw(:,no2_vis.no2_dscd_err)./1e15;

        out_array(:,10)=qdoas_raw(:,no2_vis.o4_dscd)./1e40;
        out_array(:,11)=qdoas_raw(:,no2_vis.o4_dscd_err)./1e40;

        out_array(:,12)=qdoas_raw(:,no2_vis.no2a_dscd)./1e15;
        out_array(:,13)=qdoas_raw(:,no2_vis.no2a_dscd_err)./1e15;

        out_array(:,14)=qdoas_raw(:,no2_vis.o3_dscd)./1e20;
        out_array(:,15)=qdoas_raw(:,no2_vis.o3_dscd_err)./1e20;

        out_array(:,16)=qdoas_raw(:,no2_vis.h2o_dscd)./1e23;
        out_array(:,17)=qdoas_raw(:,no2_vis.h2o_dscd_err)./1e23;

        out_array(:,18)=qdoas_raw(:,no2_vis.ring);
        out_array(:,19)=qdoas_raw(:,no2_vis.ring_err);

        % fit parameters and intensities
        out_array(:,20)=qdoas_raw(:,no2_vis.rms);
        out_array(:,21)=qdoas_raw(:,no2_vis.shift);
        out_array(:,22)=ri_440;
        out_array(:,23)=ci_1;
        out_array(:,24)=qdoas_raw(:,no2_vis.offset_0);

        %%% Write header 

        % open output file
        fid = fopen(f_out, 'w');

        % write necessary info (specific to retrieved trace gas)
        fprintf(fid(1), '%s\n', '* NofHeaderlines: 44');
        fprintf(fid(1), '%s\n', '* NofColumns: 24 (if any info missing, put -999, even if it is the whole column) ');
        fprintf(fid(1), '%s\n', '* Instrument identifier: UToronto_MAXDOAS');
        fprintf(fid(1), '%s\n', '* Retrieval code: QDOAS (v2.111, April 2016)');
        fprintf(fid(1), '%s\n', '* Created by: Kristof Bognar');
        fprintf(fid(1), '%s\n', '* Version: NO2vis_v3');
        fprintf(fid(1), '%s\n', '* X-Axis (Col 1) = Day of year (DOY) 2016 (please start with 0.0 for January 1st, 0:00 UTC)');
        fprintf(fid(1), '%s\n', '* Y1-Axis (Col 2) = Time of day in hours (UTC)');
        fprintf(fid(1), '%s\n', '* Y2-Axis (Col 3) = Total Integration Time(s) ');
        fprintf(fid(1), '%s\n', '* Y3-Axis (Col 4) = Solar Zenith Angle (°)');
        fprintf(fid(1), '%s\n', '* Y4-Axis (Col 5) = Solar Azimuth Angle (°) North=0, East=90');
        fprintf(fid(1), '%s\n', '* Y5-Axis (Col 6) = Elevation Angle (°)');
        fprintf(fid(1), '%s\n', '* Y6-Axis (Col 7) = Viewing Angle (°) North=0, East=90');

        % ... all QDOAS products for given retrieval (species+error)
        fprintf(fid(1), '%s\n', '* Y7-Axis (Col  8) = NO2_DSCD_298 (1*10^15 molec/cm2)');
        fprintf(fid(1), '%s\n', '* Y8-Axis (Col  9) = NO2_DSCD_298_Error (1*10^15 molec/cm2)');
        fprintf(fid(1), '%s\n', '* Y9-Axis (Col 10) = O4_DSCD (1*10^40 molec2/cm5)');
        fprintf(fid(1), '%s\n', '* Y10-Axis (Col 11) = O4_DSCD_Error (1*10^40 molec2/cm5)');
        fprintf(fid(1), '%s\n', '* Y11-Axis (Col 12) = NO2a_DSCD_220  (1*10^15 molec/cm2)  (Fit results for the "cold NO2 residue")');
        fprintf(fid(1), '%s\n', '* Y12-Axis (Col 13) = NO2a_DSCD_220_Error  (1*10^15 molec/cm2)');
        fprintf(fid(1), '%s\n', '* Y13-Axis (Col 14) = O3_DSCD_223 (1*10^20 molecules/cm2)');
        fprintf(fid(1), '%s\n', '* Y14-Axis (Col 15) = O3_DSCD_223_Error (1*10^20 molecules/cm2)');
        fprintf(fid(1), '%s\n', '* Y15-Axis (Col 16) = H2O_DSCD (1*10^23 molec/cm2)');
        fprintf(fid(1), '%s\n', '* Y16-Axis (Col 17) = H2O_DSCD_Error (1*10^23 molec/cm2)');
        fprintf(fid(1), '%s\n', '* Y17-Axis (Col 18) = Ring');
        fprintf(fid(1), '%s\n', '* Y18-Axis (Col 19) = Ring_Error');
        fprintf(fid(1), '%s\n', '* Y19-Axis (Col 20) = Fit RMS (in OD)');
        fprintf(fid(1), '%s\n', '* Y20-Axis (Col 21) = Spectrum shift (nm, against FRS reference)');
        fprintf(fid(1), '%s\n', '* Y21-Axis (Col 22) = Relative Intensity (counts/integration time @ 440nm)');
        fprintf(fid(1), '%s\n', '* Y22-Axis (Col 23) = Colour index: (425 nm / 440 nm)');
        fprintf(fid(1), '%s\n', '* Y23-Axis (Col 24) = Intensity offset with normalisation by I, I is the mean intensity in the spectral analysis windows, constant term');

        % other details
        fprintf(fid(1), '%s\n', '* Fit settings: 1');
        fprintf(fid(1), '%s\n', '* Fitting Window: 425-490 nm');
        fprintf(fid(1), '%s\n', '* Polynomial: 5 (6 coefficients)');
        fprintf(fid(1), '%s\n', '* Offset: Zeroth order (constant)');
        fprintf(fid(1), '%s\n', '* Calibration: Based on reference SAO solar spectra (Chance and Kurucz, 2010) --> sao2010_solref_air.dat');
        fprintf(fid(1), '%s\n', '* Wavelength adjustment: all spectra shifted and stretched against reference spectrum');
        fprintf(fid(1), '%s\n', '* Reference: zenith-sky sequential reference interpolated to the time of each individual measurement');

        % ... all the cross-section files go here
        fprintf(fid(1), '%s\n', '* NO2_298  : Vandaele et al. (1998), 298 K with I0 correction (1*10^15 molecules/cm2)  --> file: no2_298K_vanDaele.xs');
        fprintf(fid(1), '%s\n', '*  NO2_220    :  Vandaele  et  al.  (1998),  220  K  with  I0  correction  (1*10^15  molecules/cm2)  pre-orthogonalized --> file: no2a_220p298K_vanDaele_425-490.xs');
        fprintf(fid(1), '%s\n', '*   O3     : Serdyuchenko  et al., (2014), 223 K with I0 correction (1*10^20 molecules/cm2) --> file: o3_223K_SDY_air.xs ');
        fprintf(fid(1), '%s\n', '*   O4     : Thalman and Volkamer 2013, 293 K  --> file: o4_thalman_volkamer_293K_air_corrected.xs ');
        fprintf(fid(1), '%s\n', '*   H2O: HITEMP, Rothman et al., 2010 --> file: H2O_HITEMP_2010_390-700_296K_1013mbar_air.xs');
        fprintf(fid(1), '%s\n', '*      RING      : High  Resolution  calculation  with  QDOAS  according  to  Chance  and  Spurr  (1997)  and normalized as in Wagner et al. (2009) --> file: Ring_QDOAScalc_HighResSAO2010_Norm.xs');

        % ... write column headers here
        head=sprintf('*DOY\tUTC\tTint\tSZA\tSAA\tElev\tViewing_angle\tNO2_DSCD_298\tNO2_DSCD_298_error\tO4_DSCD\tO4_DSCD_error\tNO2a_DSCD_220\tNO2a_DSCD_220_Error\tO3_DSCD_223\tO3_DSCD_223_Error\tH2O_DSCD\tH2O_DSCD_Error\tRing\tRing_Error\tRMS\tSpectrum_shift\tIntens(440)\tCI(425/440)\toffset_cst');
        fprintf(fid(1), '%s\n', head);

        fclose(fid);

        %%% Writre data in the file

        dlmwrite(f_out,out_array,'delimiter','\t','precision',8,'-append')

    
    elseif selection == 2 % NO2 and O4 alt
        
        % set array size and output file name
        colnum=24;
        f_out=['CINDI_files_scanref/UToronto_MAXDOAS_36_NO2visSmall_CINDI2_' date '_v3.asc'];
        
        % create output array
        out_array=zeros(size(qdoas_raw,1),colnum);
        
         % Assign columns: general info -- same for all species
%         out_array(:,1)=day;
        out_array(:,1)=qdoas_raw(:,no2_vis.fd)-1;        
        out_array(:,2)=qdoas_raw(:,no2_vis.ft);
        out_array(:,3)=tint;
        out_array(:,4)=qdoas_raw(:,no2_vis.sza);
        out_array(:,5)=qdoas_raw(:,no2_vis.saa)+180.0;
        out_array(:,6)=qdoas_raw(:,no2_vis.elev);
        out_array(:,7)=qdoas_raw(:,no2_vis.az);

        % DSCDs and errors
        out_array(:,8)=qdoas_raw(:,no2_alt.no2_dscd)./1e15;
        out_array(:,9)=qdoas_raw(:,no2_alt.no2_dscd_err)./1e15;

        out_array(:,10)=qdoas_raw(:,no2_alt.o4_dscd)./1e40;
        out_array(:,11)=qdoas_raw(:,no2_alt.o4_dscd_err)./1e40;

        out_array(:,12)=qdoas_raw(:,no2_alt.no2a_dscd)./1e15;
        out_array(:,13)=qdoas_raw(:,no2_alt.no2a_dscd_err)./1e15;

        out_array(:,14)=qdoas_raw(:,no2_alt.o3_dscd)./1e20;
        out_array(:,15)=qdoas_raw(:,no2_alt.o3_dscd_err)./1e20;

        out_array(:,16)=qdoas_raw(:,no2_alt.h2o_dscd)./1e23;
        out_array(:,17)=qdoas_raw(:,no2_alt.h2o_dscd_err)./1e23;

        out_array(:,18)=qdoas_raw(:,no2_alt.ring);
        out_array(:,19)=qdoas_raw(:,no2_alt.ring_err);

        % fit parameters and intensities
        out_array(:,20)=qdoas_raw(:,no2_alt.rms);
        out_array(:,21)=qdoas_raw(:,no2_alt.shift);
        out_array(:,22)=ri_440;
        out_array(:,23)=ci_2;
        out_array(:,24)=qdoas_raw(:,no2_alt.offset_0);

        %%% Write header 

        % open output file
        fid = fopen(f_out, 'w');

        % write necessary info (specific to retrieved trace gas)
        fprintf(fid(1), '%s\n', '* NofHeaderlines: 44');
        fprintf(fid(1), '%s\n', '* NofColumns: 24 (if any info missing, put -999, even if it is the whole column) ');
        fprintf(fid(1), '%s\n', '* Instrument identifier: UToronto_MAXDOAS');
        fprintf(fid(1), '%s\n', '* Retrieval code: QDOAS (v2.111, April 2016)');
        fprintf(fid(1), '%s\n', '* Created by: Kristof Bognar');
        fprintf(fid(1), '%s\n', '* Version: NO2visSmall_v3');
        fprintf(fid(1), '%s\n', '* X-Axis (Col 1) = Day of year (DOY) 2016 (please start with 0.0 for January 1st, 0:00 UTC)');
        fprintf(fid(1), '%s\n', '* Y1-Axis (Col 2) = Time of day in hours (UTC)');
        fprintf(fid(1), '%s\n', '* Y2-Axis (Col 3) = Total Integration Time(s) ');
        fprintf(fid(1), '%s\n', '* Y3-Axis (Col 4) = Solar Zenith Angle (°)');
        fprintf(fid(1), '%s\n', '* Y4-Axis (Col 5) = Solar Azimuth Angle (°) North=0, East=90');
        fprintf(fid(1), '%s\n', '* Y5-Axis (Col 6) = Elevation Angle (°)');
        fprintf(fid(1), '%s\n', '* Y6-Axis (Col 7) = Viewing Angle (°) North=0, East=90');

        % ... all QDOAS products for given retrieval (species+error)
        fprintf(fid(1), '%s\n', '* Y7-Axis (Col  8) = NO2_DSCD_298 (1*10^15 molec/cm2)');
        fprintf(fid(1), '%s\n', '* Y8-Axis (Col  9) = NO2_DSCD_298_Error (1*10^15 molec/cm2)');
        fprintf(fid(1), '%s\n', '* Y9-Axis (Col 10) = O4_DSCD (1*10^40 molec2/cm5)');
        fprintf(fid(1), '%s\n', '* Y10-Axis (Col 11) = O4_DSCD_Error (1*10^40 molec2/cm5)');
        fprintf(fid(1), '%s\n', '* Y11-Axis (Col 12) = NO2a_DSCD_220  (1*10^15 molec/cm2)  (Fit results for the "cold NO2 residue")');
        fprintf(fid(1), '%s\n', '* Y12-Axis (Col 13) = NO2a_DSCD_220_Error  (1*10^15 molec/cm2)');
        fprintf(fid(1), '%s\n', '* Y13-Axis (Col 14) = O3_DSCD_223 (1*10^20 molecules/cm2)');
        fprintf(fid(1), '%s\n', '* Y14-Axis (Col 15) = O3_DSCD_223_Error (1*10^20 molecules/cm2)');
        fprintf(fid(1), '%s\n', '* Y15-Axis (Col 16) = H2O_DSCD (1*10^23 molec/cm2)');
        fprintf(fid(1), '%s\n', '* Y16-Axis (Col 17) = H2O_DSCD_Error (1*10^23 molec/cm2)');
        fprintf(fid(1), '%s\n', '* Y17-Axis (Col 18) = Ring');
        fprintf(fid(1), '%s\n', '* Y18-Axis (Col 19) = Ring_Error');
        fprintf(fid(1), '%s\n', '* Y19-Axis (Col 20) = Fit RMS (in OD)');
        fprintf(fid(1), '%s\n', '* Y20-Axis (Col 21) = Spectrum shift (nm, against FRS reference)');
        fprintf(fid(1), '%s\n', '* Y21-Axis (Col 22) = Relative Intensity (counts/integration time @ 440nm)');
        fprintf(fid(1), '%s\n', '* Y22-Axis (Col 23) = Colour index: (412 nm / 440 nm)');
        fprintf(fid(1), '%s\n', '* Y23-Axis (Col 24) = Intensity offset with normalisation by I, I is the mean intensity in the spectral analysis windows, constant term');

        % other details
        fprintf(fid(1), '%s\n', '* Fit settings: 1');
        fprintf(fid(1), '%s\n', '* Fitting Window: 411-445 nm');
        fprintf(fid(1), '%s\n', '* Polynomial: 4 (5 coefficients)');
        fprintf(fid(1), '%s\n', '* Offset: Zeroth order (constant)');
        fprintf(fid(1), '%s\n', '* Calibration: Based on reference SAO solar spectra (Chance and Kurucz, 2010) --> sao2010_solref_air.dat');
        fprintf(fid(1), '%s\n', '* Wavelength adjustment: all spectra shifted and stretched against reference spectrum');
        fprintf(fid(1), '%s\n', '* Reference: zenith-sky sequential reference interpolated to the time of each individual measurement');

        % ... all the cross-section files go here
        fprintf(fid(1), '%s\n', '* NO2_298  : Vandaele et al. (1998), 298 K with I0 correction (1*10^15 molecules/cm2)  --> file: no2_298K_vanDaele.xs');
        fprintf(fid(1), '%s\n', '*  NO2_220    :  Vandaele  et  al.  (1998),  220  K  with  I0  correction  (1*10^15  molecules/cm2)  pre-orthogonalized --> file: no2a_220p298K_vanDaele_410-446nm.xs');
        fprintf(fid(1), '%s\n', '*   O3     : Serdyuchenko  et al., (2014), 223 K with I0 correction (1*10^20 molecules/cm2) --> file: o3_223K_SDY_air.xs ');
        fprintf(fid(1), '%s\n', '*   O4     : Thalman and Volkamer 2013, 293 K  --> file: o4_thalman_volkamer_293K_air_corrected.xs ');
        fprintf(fid(1), '%s\n', '*   H2O: HITEMP, Rothman et al., 2010 --> file: H2O_HITEMP_2010_390-700_296K_1013mbar_air.xs');
        fprintf(fid(1), '%s\n', '*      RING      : High  Resolution  calculation  with  QDOAS  according  to  Chance  and  Spurr  (1997)  and normalized as in Wagner et al. (2009) --> file: Ring_QDOAScalc_HighResSAO2010_Norm.xs');

        % ... write column headers here
        head=sprintf('*DOY\tUTC\tTint\tSZA\tSAA\tElev\tViewing_angle\tNO2_DSCD_298\tNO2_DSCD_298_error\tO4_DSCD\tO4_DSCD_error\tNO2a_DSCD_220\tNO2a_DSCD_220_Error\tO3_DSCD_223\tO3_DSCD_223_Error\tH2O_DSCD\tH2O_DSCD_Error\tRing\tRing_Error\tRMS\tSpectrum_shift\tIntens(440)\tCI(412/440)\toffset_cst');
        fprintf(fid(1), '%s\n', head);

        fclose(fid);

        %%% Writre data in the file

        dlmwrite(f_out,out_array,'delimiter','\t','precision',8,'-append')       
    end
end

end
