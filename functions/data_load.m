%% Load all vintages

DATA = dir(fullfile(pathdata,'/vintages_2007_2010.xlsx'));
[status,sheets] = xlsfinfo(DATA.name);
sheets = sheetnames(DATA.name);
length_sheets_vintages_2007_2010=length(sheets);

DATA = dir(fullfile(pathdata,'/vintages_2011_2014.xlsx'));
[status,sheets] = xlsfinfo(DATA.name);
sheets = sheetnames(DATA.name);
length_sheets_vintages_2011_2014=length(sheets);

DATA = dir(fullfile(pathdata,'/vintages_2015_2018.xlsx'));
[status,sheets] = xlsfinfo(DATA.name);
sheets = sheetnames(DATA.name);
length_sheets_vintages_2015_2018=length(sheets);

DATA = dir(fullfile(pathdata,'/vintages_2019_2022.xlsx'));
[status,sheets] = xlsfinfo(DATA.name);
sheets = sheetnames(DATA.name);
length_sheets_vintages_2019_2022=length(sheets);

fs = matlab.io.datastore.FileSet("vintages_2007_2010.xlsx");
vintages_2007_2010 = spreadsheetDatastore(fs);
fs = matlab.io.datastore.FileSet("vintages_2011_2014.xlsx");
vintages_2011_2014 = spreadsheetDatastore(fs);
fs = matlab.io.datastore.FileSet("vintages_2015_2018.xlsx");
vintages_2015_2018 = spreadsheetDatastore(fs);
fs = matlab.io.datastore.FileSet("vintages_2019_2022.xlsx");
vintages_2019_2022 = spreadsheetDatastore(fs);

no_of_vintages = length_sheets_vintages_2007_2010+length_sheets_vintages_2011_2014 + length_sheets_vintages_2015_2018 + length_sheets_vintages_2019_2022;
Y_raw_vintages=cell(no_of_vintages,1);

vintage = vintages_2007_2010;
for vint_i = 1:length_sheets_vintages_2007_2010
    vintage.Sheets=vint_i;
    Y_raw_vintages{vint_i,1} = read(vintage);
end
vintage = vintages_2011_2014;
for vint_i = 1:length_sheets_vintages_2011_2014
    vintage.Sheets=vint_i;
    Y_raw_vintages{length_sheets_vintages_2011_2014 + vint_i,1} = read(vintage);
end
vintage = vintages_2015_2018;
for vint_i = 1:length_sheets_vintages_2015_2018
    vintage.Sheets=vint_i;
    Y_raw_vintages{length_sheets_vintages_2007_2010 + length_sheets_vintages_2011_2014 + vint_i,1} = read(vintage);
end
vintage = vintages_2019_2022;
for vint_i = 1:length_sheets_vintages_2019_2022
    vintage.Sheets=vint_i;
    Y_raw_vintages{length_sheets_vintages_2007_2010 + length_sheets_vintages_2011_2014 + length_sheets_vintages_2015_2018 + vint_i,1} = read(vintage);
end

% fs = matlab.io.datastore.FileSet("spf_fe_dataset.xlsx");
% spf_dataset = spreadsheetDatastore(fs);
% spf_dataset.Sheets=1;
% spf_dataset_h1 = read(spf_dataset);
% spf_dataset.Sheets=2;
% spf_dataset_h4 = read(spf_dataset);
% spf_dataset.Sheets=3;
% spf_dataset_SSP = read(spf_dataset);
% vintage.Sheets=length_sheets_vintages_2015_2022;
% Yraw = read(vintage);
% for vint_i in 1:size(Y_raw_vintages,1)


