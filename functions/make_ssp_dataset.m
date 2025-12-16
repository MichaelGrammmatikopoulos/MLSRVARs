% Make SSP dataset

% Specify the directory where your Excel files are located
srcdir = append(cd_path,'/data/spf_data/');

% Get a list of all Excel files in the directory
srcfiles = dir(fullfile(srcdir, '*.xlsx'));

% Initialize an empty structure to store the data
data_ssp = struct();

for k = 1:numel(srcfiles)
    % Get the full path to the file
    filename = fullfile(srcdir, srcfiles(k).name);
    
    % Get the names of all sheets in the file
    sheets = sheetnames(filename);

    % Convert the filename to a valid field name
    validFileName = matlab.lang.makeValidName(srcfiles(k).name);
    
    % Loop over each sheet
    for s = 1:numel(sheets)
        % Read the data from the sheet into a table
        T = readtable(filename, 'Sheet', sheets(s));
        
        % Convert the table to a structure and store it in the 'data' structure
        data_ssp.(validFileName).(sheets(s)) = table2struct(T, 'ToScalar', true);
    end
end

dates_id = [data_ssp.meanGrowth_xlsx.CPROF.YEAR data_ssp.meanGrowth_xlsx.CPROF.QUARTER];
dates_id2 = (dates_id(:,1)>=min_year & dates_id(:,1)<=2019);

% spreads_table = table([median([table2array(Yraw_table(end-120:end,5))-table2array(Yraw_table(end-120:end,15))])
% median([table2array(Yraw_table(end-120:end,5))-table2array(Yraw_table(end-120:end,16))])
% median([table2array(Yraw_table(end-120:end,5))-table2array(Yraw_table(end-120:end,17))])
% median([table2array(Yraw_table(end-120:end,5))-table2array(Yraw_table(end-120:end,18))])
% median([table2array(Yraw_table(end-120:end,5))-table2array(Yraw_table(end-120:end,19))])
% median([table2array(Yraw_table(end-120:end,5))-table2array(Yraw_table(end-120:end,20))])
% median([table2array(Yraw_table(end-120:end,5))-table2array(Yraw_table(end-120:end,21))])]);
% spreads_table(:,2)=append((Yraw_table.Properties.VariableNames(15:21)),'-FED')';
% spreads_table.Properties.VariableNames = {'Value','Spread'};
% 
spf_dataset_SSP = struct();
spf_dataset_SSP.CPIAUCSL = 100*((1+data_ssp.meanLevel_xlsx.CPI.CPI6/100).^(1/4)-1);
spf_dataset_SSP.GDP = 100*((1+data_ssp.meanGrowth_xlsx.RGDP.drgdp6/100).^(1/4)-1);
spf_dataset_SSP.UNRATE = data_ssp.meanLevel_xlsx.UNEMP.UNEMP6;
spf_dataset_SSP.PPIACO = 100*((1+data_ssp.meanLevel_xlsx.CPI.CPI6./100).^(1/4)-1);
spf_dataset_SSP.PINCOME = 100*((1+data_ssp.meanGrowth_xlsx.NGDP.dngdp6./100).^(1/4)-1);
spf_dataset_SSP.PAYEMS = 100*((1+data_ssp.meanGrowth_xlsx.EMP_PCG.dempb7./100).^(1/4)-1);
spf_dataset_SSP.INDPRO = 100*((1+data_ssp.meanGrowth_xlsx.INDPROD.dindprod6./100).^(1/4)-1);
spf_dataset_SSP.HOUST = log(data_ssp.meanLevel_xlsx.HOUSING.HOUSING6*1000);
spf_dataset_SSP.PCECC96 = 100*((1+data_ssp.meanGrowth_xlsx.RCONSUM.drconsum6/100).^(1/4)-1);
spf_dataset_SSP.M2SL = 0*spf_dataset_SSP.PCECC96;
spf_dataset_SSP.MCUMFN = 79.54*spf_dataset_SSP.PCECC96./spf_dataset_SSP.PCECC96;
spf_dataset_SSP.EXUSUK = 0*spf_dataset_SSP.PCECC96;
spf_dataset_SSP.FEDFUNDS = data_ssp.meanLevel_xlsx.TBILL.TBILL6;
spf_dataset_SSP.TB3MS = data_ssp.meanLevel_xlsx.TBILL.TBILL6;
spf_dataset_SSP.TB6MS = spf_dataset_SSP.TB3MS + 0.2;
spf_dataset_SSP.GS1   = spf_dataset_SSP.TB6MS + 0.2;
spf_dataset_SSP.GS3   = spf_dataset_SSP.GS1 + 0.2;
spf_dataset_SSP.GS5   = spf_dataset_SSP.GS3 + 0.2;
spf_dataset_SSP.GS10 = data_ssp.meanLevel_xlsx.TBOND.TBOND6;
spf_dataset_SSP.BAA = data_ssp.meanLevel_xlsx.BAABOND.BAABOND7;

% spf_dataset_SSP.CPIAUCSL = 2.5*data_ssp.meanLevel_xlsx.CPI.CPI6./data_ssp.meanLevel_xlsx.CPI.CPI6/4;
% spf_dataset_SSP.GDP = spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.UNRATE = 5*spf_dataset_SSP.CPIAUCSL./spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.PPIACO = spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.PINCOME = 4*spf_dataset_SSP.CPIAUCSL./spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.PAYEMS = spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.INDPRO = spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.HOUST = 7*spf_dataset_SSP.CPIAUCSL./spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.PCECC96 = spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.M2SL = 0*spf_dataset_SSP.CPIAUCSL./spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.MCUMFN = 79.*spf_dataset_SSP.CPIAUCSL./spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.EXUSUK = 0*spf_dataset_SSP.CPIAUCSL./spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.FEDFUNDS = 2.5*spf_dataset_SSP.CPIAUCSL./spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.TB3MS = 2.5*spf_dataset_SSP.CPIAUCSL./spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.TB6MS = spf_dataset_SSP.TB3MS + 0.2;
% spf_dataset_SSP.GS1   = spf_dataset_SSP.TB6MS + 0.2;
% spf_dataset_SSP.GS3   = spf_dataset_SSP.GS1 + 0.2;
% spf_dataset_SSP.GS5   = spf_dataset_SSP.GS3 + 0.2;
% spf_dataset_SSP.GS10 = 5*spf_dataset_SSP.CPIAUCSL./spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.BAA = 7*spf_dataset_SSP.CPIAUCSL./spf_dataset_SSP.CPIAUCSL;

% spf_dataset_SSP.CPIAUCSL    = 100*((1+data_ssp.meanLevel_xlsx.CPI10.CPI10/100).^(1/4)-1);
% spf_dataset_SSP.GDP         = 100*((1+data_ssp.meanLevel_xlsx.RGDP10.RGDP10f/100).^(1/4)-1);
% spf_dataset_SSP.UNRATE      = data_ssp.meanLevel_xlsx.UBAR.UBARf;
% spf_dataset_SSP.PPIACO      = spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.PINCOME     = spf_dataset_SSP.GDP;
% spf_dataset_SSP.PAYEMS      = spf_dataset_SSP.GDP;
% spf_dataset_SSP.INDPRO      = spf_dataset_SSP.GDP;
% spf_dataset_SSP.HOUST       = log(data_ssp.meanLevel_xlsx.HOUSING.HOUSING6*1000);
% spf_dataset_SSP.PCECC96     = spf_dataset_SSP.GDP;
% spf_dataset_SSP.M2SL        = 0*spf_dataset_SSP.GDP;
% spf_dataset_SSP.MCUMFN      = 79.54*spf_dataset_SSP.GDP./spf_dataset_SSP.GDP;
% spf_dataset_SSP.EXUSUK      = 0*spf_dataset_SSP.GDP;
% spf_dataset_SSP.FEDFUNDS    = data_ssp.meanLevel_xlsx.BILL10.BILL10f;
% spf_dataset_SSP.TB3MS       = data_ssp.meanLevel_xlsx.BILL10.BILL10f;
% spf_dataset_SSP.TB6MS       = spf_dataset_SSP.TB3MS  + 0.2;
% spf_dataset_SSP.GS1         = spf_dataset_SSP.TB6MS  + 0.2;
% spf_dataset_SSP.GS3         = spf_dataset_SSP.GS1 + 0.2;
% spf_dataset_SSP.GS5         = spf_dataset_SSP.GS3  + 0.2;
% spf_dataset_SSP.GS10        = data_ssp.meanLevel_xlsx.BOND10.BOND10f;
% spf_dataset_SSP.BAA         = data_ssp.meanLevel_xlsx.BOND10.BOND10f+2;

% Get the field names of the struct
fields = fieldnames(spf_dataset_SSP);

% Convert the first dataset to a table
i=1;
Tbl = array2table(spf_dataset_SSP.(fields{i}), 'VariableNames', {fields{i}});

% Loop over each remaining field
for i = 2:numel(fields)
    % Convert the data for the current field to a table
    data = array2table(spf_dataset_SSP.(fields{i}), 'VariableNames', {fields{i}});
    
    % Append the data to the table
    Tbl = [Tbl, data];
end
spf_dataset_SSP_2 = struct2table(spf_dataset_SSP);
spf_dataset_SSP_2 = [array2table(data_ssp.meanGrowth_xlsx.CPROF.YEAR(dates_id2,:)) spf_dataset_SSP_2(dates_id2,:)];

spf_dataset_SSP = spf_dataset_SSP_2(1:end-1,:);

% fill NAs
spf_dataset_SSP.BAA(1:20) = spf_dataset_SSP.GS10(1:20) + 2;


% spf_dataset_SSP = struct();
% spf_dataset_SSP.CPIAUCSL = 100*((1+data_ssp.meanLevel_xlsx.CPI.CPI6/100).^(1/4)-1);
% spf_dataset_SSP.GDP = 100*((1+data_ssp.meanGrowth_xlsx.RGDP.drgdp6/100).^(1/4)-1);
% spf_dataset_SSP.UNRATE = data_ssp.meanLevel_xlsx.UNEMP.UNEMP6;
% spf_dataset_SSP.PPIACO = 100*((1+data_ssp.meanLevel_xlsx.CPI.CPI6./100).^(1/4)-1);
% spf_dataset_SSP.PINCOME = 100*((1+data_ssp.meanGrowth_xlsx.NGDP.dngdp6./100).^(1/4)-1);
% spf_dataset_SSP.PAYEMS = 100*((1+data_ssp.meanGrowth_xlsx.EMP_PCG.dempb7./100).^(1/4)-1);
% spf_dataset_SSP.INDPRO = 100*((1+data_ssp.meanGrowth_xlsx.INDPROD.dindprod6./100).^(1/4)-1);
% spf_dataset_SSP.HOUST = log(data_ssp.meanLevel_xlsx.HOUSING.HOUSING6*1000);
% spf_dataset_SSP.PCECC96 = 100*((1+data_ssp.meanGrowth_xlsx.RCONSUM.drconsum6/100).^(1/4)-1);
% spf_dataset_SSP.M2SL = 0*spf_dataset_SSP.PCECC96;
% spf_dataset_SSP.MCUMFN = 79.54*spf_dataset_SSP.PCECC96./spf_dataset_SSP.PCECC96;
% spf_dataset_SSP.EXUSUK = 0*spf_dataset_SSP.PCECC96;
% spf_dataset_SSP.FEDFUNDS = data_ssp.meanLevel_xlsx.TBILL.TBILL6;
% spf_dataset_SSP.TB3MS = data_ssp.meanLevel_xlsx.TBILL.TBILL6;
% spf_dataset_SSP.TB6MS = spf_dataset_SSP.TB3MS + 0.2;
% spf_dataset_SSP.GS1   = spf_dataset_SSP.TB6MS + 0.2;
% spf_dataset_SSP.GS3   = spf_dataset_SSP.GS1 + 0.2;
% spf_dataset_SSP.GS5   = spf_dataset_SSP.GS3 + 0.2;
% spf_dataset_SSP.GS10 = data_ssp.meanLevel_xlsx.TBOND.TBOND6;
% spf_dataset_SSP.BAA = data_ssp.meanLevel_xlsx.BAABOND.BAABOND7;

% spf_dataset_SSP = struct();
% spf_dataset_SSP.CPIAUCSL    = 100*((1+data_ssp.meanLevel_xlsx.CPI10.CPI10/100).^(1/4)-1);
% spf_dataset_SSP.GDP         = 100*((1+data_ssp.meanLevel_xlsx.RGDP10.RGDP10f/100).^(1/4)-1);
% spf_dataset_SSP.UNRATE      = data_ssp.meanLevel_xlsx.UBAR.UBARf;
% spf_dataset_SSP.PPIACO      = spf_dataset_SSP.CPIAUCSL;
% spf_dataset_SSP.PINCOME     = 100*((1+data_ssp.meanGrowth_xlsx.NGDP.dngdp6./100).^(1/4)-1);
% spf_dataset_SSP.PAYEMS      = 100*((1+data_ssp.meanGrowth_xlsx.EMP_PCG.dempb7./100).^(1/4)-1);
% spf_dataset_SSP.INDPRO      = 100*((1+data_ssp.meanLevel_xlsx.RGDP10.RGDP10f/100).^(1/4)-1);
% spf_dataset_SSP.HOUST       = log(data_ssp.meanLevel_xlsx.HOUSING.HOUSING6*1000);
% spf_dataset_SSP.PCECC96     = 100*((1+data_ssp.meanGrowth_xlsx.RCONSUM.drconsum6/100).^(1/4)-1);
% spf_dataset_SSP.M2SL        = 0*spf_dataset_SSP.PCECC96;
% spf_dataset_SSP.MCUMFN      = 79.54*spf_dataset_SSP.PCECC96./spf_dataset_SSP.PCECC96;
% spf_dataset_SSP.EXUSUK      = 0*spf_dataset_SSP.PCECC96;
% spf_dataset_SSP.FEDFUNDS    = data_ssp.meanLevel_xlsx.BILL10.BILL10f;
% spf_dataset_SSP.TB3MS       = data_ssp.meanLevel_xlsx.BILL10.BILL10f;
% spf_dataset_SSP.TB6MS       = spf_dataset_SSP.TB3MS  + 0.15;
% spf_dataset_SSP.GS1         = spf_dataset_SSP.TB6MS  + 0.15;
% spf_dataset_SSP.GS3         = spf_dataset_SSP.GS1 + 0.15;
% spf_dataset_SSP.GS5         = spf_dataset_SSP.GS3  + 0.15;
% spf_dataset_SSP.GS10        = data_ssp.meanLevel_xlsx.BOND10.BOND10f;
% spf_dataset_SSP.BAA         = data_ssp.meanLevel_xlsx.BOND10.BOND10f+2;
% spf_dataset_SSP.PREMIA      = spf_dataset_SSP.GS10-0.5*(spf_dataset_SSP.GS1+spf_dataset_SSP.GS3);
