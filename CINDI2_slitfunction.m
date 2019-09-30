% generate slit function for CINDI2 spectra submission

% slf required is a single peak at any wavelength (preferrably in middle of detector)
% area of peak should be normalized to 1
% wavelength should be given as +- from peak

%% Hg lamp test spectrum from CINDI2
spec_raw = get_spectrum_v2('06091611.csv','11:03:00');

% integration time in seconds
% exp time = 14 ms
% scans = 100
tint=1.4;

%% Read and subtract DC and bias files
load dark_c.inp
load bias.inp

spec = spec_raw - dark_c' * tint - bias';

%% Get desired peak and normalize area to 1

xpk=961;
range=23;
x=xpk-range:xpk+range;

% convert to nm for 600 grating
lambda=(x-xpk).*0.107;

spec_norm = spec(x)./sum(spec(x));

% plot(lambda,spec_norm)

%% print results to file

fid=fopen('UToronto_MAXDOAS_36_channel_sltfct_CINDI2_20160924_v1.asc','w');

% peak is at 445.93 in day 253 calibration file
fprintf(fid, '%s\n', '# SLIT FUNCTION AT 446 nm');
fprintf(fid, '%s\n', '# Institute = University of Toronto');
fprintf(fid, '%s\n', '# PI name = Kimberly Strong (strong@atmosp.physics.utoronto.ca)');
fprintf(fid, '%s\n', '# Instrument = VIS');

out = [lambda; spec_norm];

fprintf(fid, '%.6f %.6f \n', out);

fclose(fid);
