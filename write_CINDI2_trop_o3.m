%%% Program to write CINDI-2 data files for resubmission using scan
%%% references -- only works if date and time fields are not in the QDOAS
%%% file

cd /home/kristof/work/GBS/PEARL-GBS/CINDI-2_2016/QDOAS_output/


% load qdoas table (imported manually)
load('trop_o3.mat');

for day=256:272
    %% Read onde day of QDOAS data
    
    ind_day=find(floor(trop_o3.Fractionalday)==day);
    
    if isempty(ind_day), continue, end
    
    % filter non-MAXDOAS data for d256 (some of the horizon scan remains,
    % but that's ok)
    if day==256
        
        ind_good=find(trop_o3.Elevviewingangle==1 | ...
                      trop_o3.Elevviewingangle==2 | ...
                      trop_o3.Elevviewingangle==3 | ...
                      trop_o3.Elevviewingangle==4 | ...
                      trop_o3.Elevviewingangle==5 | ...
                      trop_o3.Elevviewingangle==6 | ...
                      trop_o3.Elevviewingangle==8 | ...
                      trop_o3.Elevviewingangle==15 | ...
                      trop_o3.Elevviewingangle==30 ...
                      );
        
        ind_day=intersect(ind_day,ind_good);
    end
    

    % get corresponding date
    [dd,mm] = Julian2Date(2016,day);
    date=datetime(2016,mm,dd,'Format','yyyyMMdd');
    date=datestr(date,'yyyymmdd');

    disp(['***   Writing data for ' num2str(mm) '/' num2str(dd) '   ***'])

    trop_o3_tmp=trop_o3(ind_day,:);
    

    %% Get t_int and intensities

    % calculate total integration time (tint) in seconds (tint in QDOAS file is
    % seconds and it is for single spectrum)
    tint=trop_o3_tmp.Scans.*trop_o3_tmp.Tint;

    % calculate relative intensity (ri) at given lambda (use exposure for
    % single measurement here!)
    ri_500=trop_o3_tmp.Fluxes500./trop_o3_tmp.Tint;

    % calculate color index (ci, as defined in data protocol)
    ci_3=trop_o3_tmp.Fluxes440./trop_o3_tmp.Fluxes500;

    %% Assign selection-specific values

    % set array size and output file name
    colnum=27;
    f_out=['CINDI_files_trop_o3/UToronto_MAXDOAS_36_O3vis_CINDI2_seq_' date '_v1.asc'];

    % create output array
    out_array=zeros(size(trop_o3_tmp,1),colnum);

    % Assign columns: general info -- same for all species
    %         out_array(:,1)=day;
    out_array(:,1)=trop_o3_tmp.Fractionalday-1;
    out_array(:,2)=trop_o3_tmp.Fractionaltime;
    out_array(:,3)=tint;
    out_array(:,4)=trop_o3_tmp.SZA;
    out_array(:,5)=trop_o3_tmp.SolarAzimuthAngle+180.0;
    out_array(:,6)=trop_o3_tmp.Elevviewingangle;
    out_array(:,7)=trop_o3_tmp.Azimviewingangle;

    % DSCDs and errors, divided by factors as specified in header text
    out_array(:,8)=trop_o3_tmp.O3SlColo3./1e18;
    out_array(:,9)=trop_o3_tmp.O3SlErro3./1e18;

    out_array(:,10)=trop_o3_tmp.O3SlColo3_293K./1e18;
    out_array(:,11)=trop_o3_tmp.O3SlErro3_293K./1e18;

    out_array(:,12)=trop_o3_tmp.O3SlColo4./1e40;
    out_array(:,13)=trop_o3_tmp.O3SlErro4./1e40;

    out_array(:,14)=trop_o3_tmp.O3SlColno2a_220p298K./1e15;
    out_array(:,15)=trop_o3_tmp.O3SlErrno2a_220p298K./1e15;

    out_array(:,16)=trop_o3_tmp.O3SlColno2./1e15;
    out_array(:,17)=trop_o3_tmp.O3SlErrno2./1e15;

    out_array(:,18)=trop_o3_tmp.O3SlColh2o./1e23;
    out_array(:,19)=trop_o3_tmp.O3SlErrh2o./1e23;

    out_array(:,20)=trop_o3_tmp.O3SlColRing;
    out_array(:,21)=trop_o3_tmp.O3SlErrRing;


    % fit parameters and intensities
    out_array(:,22)=trop_o3_tmp.O3RMS;
    out_array(:,23)=trop_o3_tmp.O3ShiftSpectrum;

    out_array(:,24)=ri_500;
    out_array(:,25)=ci_3;

    out_array(:,26)=trop_o3_tmp.O3OffsetConstant;
    out_array(:,27)=trop_o3_tmp.O3OffsetOrder1;

    %%% Write header 

    % open output file
    fid = fopen(f_out, 'w');

    % write necessary info (specific to retrieved trace gas)
    fprintf(fid(1), '%s\n', '* NofHeaderlines: 48');
    fprintf(fid(1), '%s\n', '* NofColumns: 27 (if any info missing, put -999, even if it is the whole column) ');
    fprintf(fid(1), '%s\n', '* Instrument identifier: UToronto_MAXDOAS');
    fprintf(fid(1), '%s\n', '* Retrieval code: QDOAS (v3.1, June 2017)');
    fprintf(fid(1), '%s\n', '* Created by: Kristof Bognar');
    fprintf(fid(1), '%s\n', '* Version: O3vis_v1');
    fprintf(fid(1), '%s\n', '* X-Axis (Col 1) = Day of year (DOY) 2016 (please start with 0.0 for January 1st, 0:00 UTC)');
    fprintf(fid(1), '%s\n', '* Y1-Axis (Col 2) = Time of day in hours (UTC)');
    fprintf(fid(1), '%s\n', '* Y2-Axis (Col 3) = Total Integration Time(s) ');
    fprintf(fid(1), '%s\n', '* Y3-Axis (Col 4) = Solar Zenith Angle (°)');
    fprintf(fid(1), '%s\n', '* Y4-Axis (Col 5) = Solar Azimuth Angle (°) North=0, East=90');
    fprintf(fid(1), '%s\n', '* Y5-Axis (Col 6) = Elevation Angle (°)');
    fprintf(fid(1), '%s\n', '* Y6-Axis (Col 7) = Viewing Angle (°) North=0, East=90');

    % ... all QDOAS products for given retrieval (species+error)
    fprintf(fid(1), '%s\n', '* Y13-Axis (Col 8) = O3_DSCD_223 (1*10^18 molecules/cm2)');
    fprintf(fid(1), '%s\n', '* Y14-Axis (Col 9) = O3_DSCD_223_Error (1*10^18 molecules/cm2)');
    fprintf(fid(1), '%s\n', '* Y13-Axis (Col 10) = O3_DSCD_293 (1*10^18 molecules/cm2)');
    fprintf(fid(1), '%s\n', '* Y14-Axis (Col 11) = O3_DSCD_293_Error (1*10^18 molecules/cm2)');
    fprintf(fid(1), '%s\n', '* Y9-Axis (Col 12) = O4_DSCD (1*10^40 molec2/cm5)');
    fprintf(fid(1), '%s\n', '* Y10-Axis (Col 13) = O4_DSCD_Error (1*10^40 molec2/cm5)');
    fprintf(fid(1), '%s\n', '* Y11-Axis (Col 14) = NO2_DSCD_220  (1*10^15 molec/cm2)  (Fit results for the "cold NO2 residue")');
    fprintf(fid(1), '%s\n', '* Y12-Axis (Col 15) = NO2_DSCD_220_Error  (1*10^15 molec/cm2)');
    fprintf(fid(1), '%s\n', '* Y7-Axis (Col  16) = NO2_DSCD_298 (1*10^15 molec/cm2)');
    fprintf(fid(1), '%s\n', '* Y8-Axis (Col  17) = NO2_DSCD_298_Error (1*10^15 molec/cm2)');
    fprintf(fid(1), '%s\n', '* Y15-Axis (Col 18) = H2O_DSCD (1*10^23 molec/cm2)');
    fprintf(fid(1), '%s\n', '* Y16-Axis (Col 19) = H2O_DSCD_Error (1*10^23 molec/cm2)');
    fprintf(fid(1), '%s\n', '* Y17-Axis (Col 20) = Ring');
    fprintf(fid(1), '%s\n', '* Y18-Axis (Col 21) = Ring_Error');
    fprintf(fid(1), '%s\n', '* Y19-Axis (Col 22) = Fit RMS (in OD)');
    fprintf(fid(1), '%s\n', '* Y20-Axis (Col 23) = Spectrum shift (nm, against FRS reference)');
    fprintf(fid(1), '%s\n', '* Y21-Axis (Col 24) = Relative Intensity (counts/integration time @ 500nm)');
    fprintf(fid(1), '%s\n', '* Y22-Axis (Col 25) = Colour index: (440 nm / 500 nm)');
    fprintf(fid(1), '%s\n', '* Y23-Axis (Col 26) = Intensity offset with normalisation by I, I is the mean intensity in the spectral analysis windows, constant term');
    fprintf(fid(1), '%s\n', '* Y23-Axis (Col 27) = Intensity offset, linear term');

    % other details
    fprintf(fid(1), '%s\n', '* Fit settings: 1');
    fprintf(fid(1), '%s\n', '* Fitting Window: 450-520 nm');
    fprintf(fid(1), '%s\n', '* Polynomial: 5 (6 coefficients)');
    fprintf(fid(1), '%s\n', '* Offset: Zeroth order (constant)');
    fprintf(fid(1), '%s\n', '* Calibration: Based on reference SAO solar spectra (Chance and Kurucz, 2010) --> sao2010_solref_air.dat');
    fprintf(fid(1), '%s\n', '* Wavelength adjustment: all spectra shifted and stretched against reference spectrum');
    fprintf(fid(1), '%s\n', '* Reference: zenith-sky sequential reference interpolated to the time of each individual measurement');

    % ... all the cross-section files go here
    fprintf(fid(1), '%s\n', '* O3_223    : Serdyuchenko  et al., (2014), 223 K with I0 correction (1*10^20 molecules/cm2) --> file: o3_223K_SDY_air.xs ');
    fprintf(fid(1), '%s\n', '* O3_293    : Serdyuchenko  et al., (2014), 293 K with I0 correction (1*10^20 molecules/cm2) --> file: o3_293K_SDY_air.xs ');
    fprintf(fid(1), '%s\n', '* NO2_298   : Vandaele et al. (1998), 298 K with I0 correction (1*10^15 molecules/cm2)  --> file: no2_298K_vanDaele.xs');
    fprintf(fid(1), '%s\n', '* NO2_220   : Vandaele  et  al.  (1998),  220  K  with  I0  correction  (1*10^15  molecules/cm2)  pre-orthogonalized --> file: no2a_220p298K_vanDaele_450-550.xs');
    fprintf(fid(1), '%s\n', '* O4        : Thalman and Volkamer 2013, 293 K  --> file: o4_thalman_volkamer_293K_air_corrected.xs ');
    fprintf(fid(1), '%s\n', '* H2O       : HITEMP, Rothman et al., 2010 --> file: H2O_HITEMP_2010_390-700_296K_1013mbar_air.xs');
    fprintf(fid(1), '%s\n', '* RING      : High  Resolution  calculation  with  QDOAS  according  to  Chance  and  Spurr  (1997)  and normalized as in Wagner et al. (2009) --> file: Ring_QDOAScalc_HighResSAO2010_Norm.xs');

    % ... write column headers here
    head=sprintf('*DOY\tUTC\tTint\tSZA\tSAA\tElev\tViewing_angle\tO3_DSCD_223\tO3_DSCD_223_Error\tO3_DSCD_293\tO3_DSCD_293_Error\tO4_DSCD\tO4_DSCD_error\tNO2_DSCD_220\tNO2_DSCD_220_Error\tNO2_DSCD_298\tNO2_DSCD_298_error\tH2O_DSCD\tH2O_DSCD_Error\tRing\tRing_Error\tRMS\tSpectrum_shift\tIntens(500)\tCI(440/500)\toffset_cst\toffset_lin');
    fprintf(fid(1), '%s\n', head);

    fclose(fid);

    %%% Writre data in the file

    dlmwrite(f_out,out_array,'delimiter','\t','precision',8,'-append')


end
