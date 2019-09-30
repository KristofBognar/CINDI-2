function CINDI_2_tg_avk_output_format(auto,version,species,data_path,data_source)
if nargin == 1
    if auto == 0
        %%%% input %%%%
        version = 'v2';
        species ='NO2_460nm';
        %%% profile data path %%%
        %data_path = 'C:\SCIATRAN2\TRACEGAS_RETRIEVAL_v-1-2\Campaign\Retrieval_settings_A\CINDI2\20160914_2\av_kernels\';
        %data_path = 'C:\SCIATRAN2\TRACEGAS_RETRIEVAL_v-1-2\Campaign\Retrieval_settings_A\20160915\av_kernels\';
        data_path = 'C:\SCIATRAN2\TRACEGAS_RETRIEVAL_v-1-2\Campaign\Retrieval_settings_A\CINDI2\no_offset\20160915\av_kernels\';
    elseif auto == 1
        disp('Warning: mising version, species, and data path info! ');
    end
elseif nargin ~= 5
    disp('Warning: mising version, species, and data path info! ');
else
    data_path = [data_path 'av_kernels/'];
end


cd(data_path);




%%%% read in data from the OEM & print%%%%
cd ..;
try
    mkdir cindi2;
end
cd(data_path);
tmp=dir();
list={tmp.name};
list(1:2)=[];
% time_stamp =list(:,5:end-4);

N = length(list);
for i = 1:1:N
    cd(data_path);
    data = importfile_avk(list{i});
    h = importfile_height(list{i});
    
    time_stamp =list{i}(5:end-4);
    
    cd ../cindi2;
    if i <= 9
        avk_name = ['U_Toronto_MAXDOAS_' data_source '_'  species '_CINDI2_' time_stamp(1:8) '_' version '_akernel_#0' num2str(i) '.asc'];
    else
        avk_name = ['U_Toronto_MAXDOAS_' data_source '_'  species '_CINDI2_' time_stamp(1:8) '_' version '_akernel_#' num2str(i) '.asc'];
    end
    dlmwrite(avk_name,data,'delimiter','\t','precision','%.5f');
    plot_avk(h,data,avk_name);
end


function plot_avk(h,data,avk_name)
save_fig = 1;
close all;
figure;hold all;
N = size(h);
for i = 1:1:N(1)
    plot(data(:,i),h);
end
xlabel(['Averaging Kernel']);
ylabel(['Altitude (km)']);
print_setting(1/4,save_fig,[avk_name(1:end-4)]);

function avk201609120608 = importfile_avk(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   AVK201609120608 = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   AVK201609120608 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data
%   from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   avk201609120608 = importfile('avk_20160912_0608.dat', 2, 22);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/11/28 10:57:39

%% Initialize variables.
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format string for each line of text:
%   column6: double (%f)
%	column7: double (%f)
%   column8: double (%f)
%	column9: double (%f)
%   column10: double (%f)
%	column11: double (%f)
%   column12: double (%f)
%	column13: double (%f)
%   column14: double (%f)
%	column15: double (%f)
%   column16: double (%f)
%	column17: double (%f)
%   column18: double (%f)
%	column19: double (%f)
%   column20: double (%f)
%	column21: double (%f)
%   column22: double (%f)
%	column23: double (%f)
%   column24: double (%f)
%	column25: double (%f)
%   column26: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*70s%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%14f%f%[^\n\r]';

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
avk201609120608 = [dataArray{1:end-1}];


function avk201609120608 = importfile_height(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   AVK201609120608 = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   AVK201609120608 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data
%   from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   avk201609120608 = importfile('avk_20160912_0608.dat', 2, 22);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/11/28 10:57:59

%% Initialize variables.
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%14f%[^\n\r]';

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
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
avk201609120608 = [dataArray{1:end-1}];



