Y_raw_Transformed_vintages         = cell(no_of_vintages,1);
adf_test_transformed_vintages      = cell(no_of_vintages,1);
quarter_names_initial              = cell(no_of_vintages,1);

for vint_i = 1:no_of_vintages

    Yraw_table = Y_raw_vintages{vint_i,1};

    %% Create inflation expectations proxy by using the spread of the 10y and 1y bonds.
    pos_LR_bond = find(ismember(Yraw_table.Properties.VariableNames,{'GS10'}));
    pos_SR_bond = find(ismember(Yraw_table.Properties.VariableNames, {'GS1'}));
    Yraw_table{:,size(Yraw_table,2)+1} = movmean(Yraw_table{:,pos_LR_bond},4)-movmean(Yraw_table{:,pos_SR_bond},4);
    Yraw_table.Properties.VariableNames{end} = 'INFEXP';
% 
%     pos_BAA_bond = find(ismember(Yraw_table.Properties.VariableNames,{'BAA'}));
%     pos_3y_bond = find(ismember(Yraw_table.Properties.VariableNames, {'GS3'}));
%     Yraw_table{:,size(Yraw_table,2)+1} = movmean(Yraw_table{:,pos_BAA_bond}-Yraw_table{:,pos_3y_bond},4);
%     Yraw_table.Properties.VariableNames{end} = 'INFEXP2';
% 
%     pos_5y_bond = find(ismember(Yraw_table.Properties.VariableNames,{'GS5'}));
%     pos_3m_bond = find(ismember(Yraw_table.Properties.VariableNames, {'TB3MS'}));
%     Yraw_table{:,size(Yraw_table,2)+1} = movmean(Yraw_table{:,pos_5y_bond}-Yraw_table{:,pos_3m_bond},4);
%     Yraw_table.Properties.VariableNames{end} = 'INFEXP3';

    % Update teh initial Y_raw
    Y_raw_vintages{vint_i,1} = Yraw_table;

    %%
    name_vint_i_last_historical = string(Y_raw_vintages{vint_i,1}.observation(end,1));

    var_mnemonic = string(Yraw_table.Properties.VariableNames);
    var_mnemonic = var_mnemonic(var_mnemonic~='observation');
   
    var_mnemonic = ["CPIAUCSL","GDP","UNRATE",...
        "PPIACO","PINCOME","PAYEMS","INDPRO","HOUST","PCECC96","M2SL","MCUMFN",'EXUSUK',...
        "FEDFUNDS","TB3MS","TB6MS","GS1","GS3","GS5","GS10","BAA","INFEXP"];
%     signs = {'<', '<', '>', '<', '<', '<', '<', '<', '<', '<', '<', '<', '>', '>', '>', '>', '>', '>' };
    ydates = Yraw_table.observation;

    Rolling_window = 0;
    if Rolling_window == 0
        Y_raw_vintages{vint_i,1} = Y_raw_vintages{vint_i,1}(:,['observation' var_mnemonic]);
        Yraw = table2array(Yraw_table(:,var_mnemonic));
    else
        Yraw = table2array(Yraw_table(vint_i:end,var_mnemonic));
    end
   
    % 7 - percent change differences
    %% Create the transformation codes
    tcode_full =cell2mat(transformation_code_creator(var_mnemonic));

    %% Transform the raw dataset to the corresponding log-diff, diff, etc.
    YrawT = Yraw;
    for jj = 1:size(YrawT,2)
        if tcode_full(1,jj)==6 % 6 - log second differences
            YrawT(:,jj) = 100*(log(Yraw(:,jj))-2*log(lagmatrix(Yraw(:,jj),1))+log(lagmatrix(Yraw(:,jj),2)));
        elseif tcode_full(1,jj)==5 % 5 - log differences
            YrawT(:,jj) = 100*(log(Yraw(:,jj))-log(lagmatrix(Yraw(:,jj),1)));
            if var_mnemonic(jj) == 'UNRATE'
               YrawT(:,jj) = YrawT(:,jj)/100;
            end
        elseif tcode_full(1,jj)==4 % 4 - log transformation
            YrawT(:,jj) = log(Yraw(:,jj));
        elseif tcode_full(1,jj)==3 % 3 - second differences
            YrawT(:,jj) = Yraw(:,jj)-2*lagmatrix(Yraw(:,jj),1)+lagmatrix(Yraw(:,jj),2);
        elseif tcode_full(1,jj)==2 % 2 - first differences
            YrawT(:,jj) = Yraw(:,jj)-lagmatrix(Yraw(:,jj),1);
        elseif tcode_full(1,jj)==1 % 1 - no transformation
            YrawT(:,jj) = Yraw(:,jj);
        end
    end

    adftest_list = [1:size(YrawT,2)];

    for j = 1:size(YrawT,2)
        adftest_list(:,j) = adftest(YrawT(:,j));
    end

    adf_test_transformed_vintages{vint_i,1} = adftest_list;
    Y_raw_Transformed_vintages{vint_i,1}=YrawT;
    quarter_names_initial{vint_i,1} = name_vint_i_last_historical;
end

Yraw_table_last_vintage=Yraw_table;

spf_dataset_SSP = [spf_dataset_SSP(:,1) spf_dataset_SSP(:,var_mnemonic)];

