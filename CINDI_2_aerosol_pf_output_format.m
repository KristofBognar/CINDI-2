function CINDI_2_aerosol_pf_output_format(auto, version, species, data_path,data_source)

if nargin == 1
    if auto == 0
        %%%% input %%%%
        version = 'v2';
        species ='O4_477nm';
        %%% profile data path %%%
        %data_path = 'C:\SCIATRAN2\AEROSOL_RETRIEVAL_v-1-2\Campaign\Retrieval_settings_A\CINDI2\20160914_2\profiles';
        %data_path = 'C:\SCIATRAN2\AEROSOL_RETRIEVAL_v-1-2\Campaign\Retrieval_settings_A\20160915\profiles\';
%         data_path ='C:\SCIATRAN2\AEROSOL_RETRIEVAL_v-1-2\Campaign\Retrieval_settings_A\CINDI2\no_offset\20160914\profiles\';
        %%% VCD %%%
        cd ..; cd general;
        file_nm_AOT = 'retrieval_20160914.dat';
        data_AOT = importfile_AOT(file_nm_AOT);
    elseif auto == 1
        disp('Warning: mising version, species, and data path info! ');
    end

elseif nargin ~= 5
    disp('Warning: mising version, species, and data path info! ');
else
    
    %%% VCD %%%
    cd ..; cd general;
    file_nm_AOT = ['retrieval_' data_path(end-8:end-1) '.dat'];
    data_AOT = importfile_AOT(file_nm_AOT);
    
    data_path = [data_path 'profiles/'];
end
            
cd(data_path);


%%%% read in data from the OEM %%%%
cd(data_path);
tmp=dir();
list={tmp.name};
list(1:2)=[];
time_stamp =list{1}(11:18);
yyyy = str2double(list{1}(11:14));
mm = str2double(list{1}(15:16));
dd = str2double(list{1}(17:18));
julian = Date2Julian(yyyy,mm,dd);
N = length(list);
for i = 1:1:N
    data = importfile_profiles(list{i});
    AOD_profiles(:,i) = data.retrieved;
    err_AOD_profiles(:,i) = data.err_retrieved;
    apriori_profiles(:,i) = data.apriori;
    HHMM = list{i}(20:23);
    fday = datenum(HHMM,'HHMM');
    fday = fday - fix(fday) + julian;
    time_profiles(:,i) = fday;
    z = data.z;
end


%%% print data to CINDI-2 format %%%
cd ..;

if ~exist('cindi2','dir'), mkdir('cindi2'), end
cd cindi2;

AOD_profile_name = ['U_Toronto_MAXDOAS_' data_source '_'  species '_CINDI2_' time_stamp '_' version '_profile.asc'];
dlmwrite(AOD_profile_name,AOD_profiles, 'delimiter','\t' ,'precision','%.8e');

err_AOD_profile_name = ['U_Toronto_MAXDOAS_' data_source '_'  species '_CINDI2_' time_stamp '_' version '_error.asc'];
dlmwrite(err_AOD_profile_name,err_AOD_profiles,'delimiter','\t','precision','%.8e');

apriori_profile_name = ['U_Toronto_MAXDOAS_' data_source '_'  species '_CINDI2_' time_stamp '_' version '_apriori_profile.asc'];
dlmwrite(apriori_profile_name,apriori_profiles,'delimiter','\t','precision','%.8e');

time_profile_name = ['U_Toronto_MAXDOAS_' data_source '_'  species '_CINDI2_' time_stamp '_' version '_time.asc'];
dlmwrite(time_profile_name,time_profiles','delimiter','\t','precision','%.8e');

altitude_profile_name = ['U_Toronto_MAXDOAS_' data_source '_'  species '_CINDI2_' time_stamp '_' version '_altitude.asc'];
altitude_profiles = [z,z+0.2,z+0.1];
dlmwrite(altitude_profile_name,altitude_profiles,'delimiter','\t','precision','%.1e');

O4_AOT_name = ['U_Toronto_MAXDOAS_' data_source '_'  species '_CINDI2_' time_stamp '_' version '_AOD.asc'];
O4_AOT = data_AOT.AOT;
dlmwrite(O4_AOT_name,O4_AOT,'delimiter','\t','precision','%.8e');




function prof477nm201609140608 = importfile_profiles(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   PROF477NM201609140608 = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   PROF477NM201609140608 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads
%   data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   prof477nm201609140608 = importfile('prof477nm_20160914_0608.dat', 2, 22);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/09/20 14:40:53

%% Initialize variables.
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%14f%14f%14f%14f%14f%14f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
prof477nm201609140608 = table(dataArray{1:end-1}, 'VariableNames', {'z','apriori','rel_err_apriori','retrieved','err_retrieved','err_smooth','err_noise'});



function retrieval20160914 = importfile_AOT(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   RETRIEVAL20160914 = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   RETRIEVAL20160914 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data
%   from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   retrieval20160914 = importfile('retrieval_20160914.dat', 2, 14);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/09/19 16:48:27

%% Initialize variables.
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: datetimes (%{dd/MM/yyyy}D)
%	column2: datetimes (%{HH:mm:ss}D)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%13{dd/MM/yyyy}D%10{HH:mm:ss}D%8f%8f%14f%14f%14f%14f%14f%14f%14f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
retrieval20160914 = table(dataArray{1:end-1}, 'VariableNames', {'date','time','Nr','n_iter','chisq','H','d_s','gamma','gamma_max','scale','AOT','err_AOT'});

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% retrieval20160914.date=datenum(retrieval20160914.date);
% retrieval20160914.time=datenum(retrieval20160914.time);

