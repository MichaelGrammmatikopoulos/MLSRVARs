%% Create tables for export to LATEX

mkdir(cd,append('results\',which_VAR,'\tables'))

table1_list_of_variables = table(pretty_names', convertStringsToChars(var_mnemonic_i'), tcode', minnesotaPriorMean);
table1_list_of_variables = renamevars(table1_list_of_variables,["Var1","Var2","Var3","minnesotaPriorMean"],["Variable","ALFRED code","Transformation code","Minnesota Prior"]);
title_ = 'List of Variables';
subtitle_ = ['Note: Data gathered from the ALFRED-MD. Includes vintages from 2008Q1 to 2019Q4. The transformations legend is the following: ' ...
    '6: log second differences, 5: log differences, 4: log transformation, 3: second differences, 2: first differences, 1: no transformation'];
table2latex(table1_list_of_variables,append('results\',which_VAR,'\tables\',which_VAR,'_','table1_list_of_variables'),title_,subtitle_,1);

format bank
% ndxMAINVARS = find(ismember(var_mnemonic, {'GDP','UNRATE','CPIAUCSL','FEDFUNDS'}));

pretty_names_short = char(var_mnemonic_i);
pretty_names_short = string((pretty_names_short(:,1:3,:)));
make_strings_for_graphs

table2_rRMSEs_h1 = string(sprintfc('%0.2f', rRMSEs_h1 ));
ndx_asterisks = rRMSEs_h1<1;
asterisks_rRMSEs_h1 = string([ndx_asterisks.*GWtests_all_models_RMSE_h1]);
asterisks_rRMSEs_h1=strrep(asterisks_rRMSEs_h1,'0','');
asterisks_rRMSEs_h1=strrep(asterisks_rRMSEs_h1,'1','*');
asterisks_rRMSEs_h1=strrep(asterisks_rRMSEs_h1,'2','**');
asterisks_rRMSEs_h1=strrep(asterisks_rRMSEs_h1,'3','***');
table2_rRMSEs_h1 = cellstr(append(table2_rRMSEs_h1,asterisks_rRMSEs_h1));
table2_rRMSEs_h1_final = array2table([cellstr(all_models_pretty) table2_rRMSEs_h1]);
table2_rRMSEs_h1_final.Properties.VariableNames = ["Model" pretty_names_short];
table2_rRMSEs_h1_final(ndxBENCH,:)=[];
title_=append(title_VAR, 'RMSE Forecast performance vs. SV-VAR, h = 1');
subtitle_=['Note: Comparison with the "SV-VAR" (Minnesota Stochastic Volatility, baseline ' ...
    'in denominator) against different specifications. Values below 1 ' ...
    'indicate improvement over baseline. Evaluation windows with forecast ' ...
    'origins from 2008Q4 through 2019:Q3. Real time vintages, from the ALFRED database. Significance is assessed by the Giacomini-White test using ' ...
    'Newey-West standard errors. Legend: SV: Stochastic Volatility, SS: Steady State, SR: Shadow Rate'];
table2latex(table2_rRMSEs_h1_final,append('results\',which_VAR,'\tables\',which_VAR,'_','table2_rRMSEs_h1_final'),title_,subtitle_,1);

table3_rRMSEs_h4 = string(sprintfc('%0.2f', rRMSEs_h4 ));
ndx_asterisks = rRMSEs_h4<1;
asterisks_rRMSEs_h4 = string([ndx_asterisks.*GWtests_all_models_RMSE_h4]);
asterisks_rRMSEs_h4=strrep(asterisks_rRMSEs_h4,'0','');
asterisks_rRMSEs_h4=strrep(asterisks_rRMSEs_h4,'1','*');
asterisks_rRMSEs_h4=strrep(asterisks_rRMSEs_h4,'2','**');
asterisks_rRMSEs_h4=strrep(asterisks_rRMSEs_h4,'3','***');
table3_rRMSEs_h4 = cellstr(append(table3_rRMSEs_h4,asterisks_rRMSEs_h4));
table3_rRMSEs_h4_final = array2table([cellstr(all_models_pretty) table3_rRMSEs_h4]);
table3_rRMSEs_h4_final.Properties.VariableNames = ["Model" pretty_names_short];
table3_rRMSEs_h4_final(ndxBENCH,:)=[];
title_=append(title_VAR,'RMSE Forecast performance vs. SV-VAR, h = 4');
subtitle_=['Note: Comparison with the "SV-VAR" (Minnesota Stochastic Volatility, baseline ' ...
    'in denominator) against different specifications. Values below 1 ' ...
    'indicate improvement over baseline. Evaluation windows with forecast ' ...
    'origins from 2008Q1 through 2018:Q4. Real time vintages, from the ALFRED database. Significance is assessed by the Giacomini-White test using ' ...
    'Newey-West standard errors. Legend: SV: Stochastic Volatility, SS: Steady State, SR: Shadow Rate'];
table2latex(table3_rRMSEs_h4_final,append('results\',which_VAR,'\tables\',which_VAR,'_','table3_rRMSEs_h4_final'),title_,subtitle_,1);

table4_rMAEs_h1 = string(sprintfc('%0.2f', rMAEs_h1 ));
ndx_asterisks = rMAEs_h1<1;
asterisks_rMAEs_h1 = string([ndx_asterisks.*GWtests_all_models_MAE_h1]);
asterisks_rMAEs_h1=strrep(asterisks_rMAEs_h1,'0','');
asterisks_rMAEs_h1=strrep(asterisks_rMAEs_h1,'1','*');
asterisks_rMAEs_h1=strrep(asterisks_rMAEs_h1,'2','**');
asterisks_rMAEs_h1=strrep(asterisks_rMAEs_h1,'3','***');
table4_rMAEs_h1 = cellstr(append(table4_rMAEs_h1,asterisks_rMAEs_h1));
table4_rMAEs_h1_final = array2table([cellstr(all_models_pretty) table4_rMAEs_h1]);
table4_rMAEs_h1_final.Properties.VariableNames = ["Model" pretty_names_short];
table4_rMAEs_h1_final(ndxBENCH,:)=[];
title_=append(title_VAR,'MAE Forecast performance vs. SV-VAR, h = 1');
subtitle_=['Note: Comparison with the "SV-VAR" (Minnesota Stochastic Volatility, baseline ' ...
    'in denominator) against different specifications. Values below 1 ' ...
    'indicate improvement over baseline. Evaluation windows with forecast ' ...
    'origins from 2008Q4 through 2019:Q3. Real time vintages, from the ALFRED database. Significance is assessed by the Giacomini-White test using ' ...
    'Newey-West standard errors. Legend: SV: Stochastic Volatility, SS: Steady State, SR: Shadow Rate'];
table2latex(table4_rMAEs_h1_final,append('results\',which_VAR,'\tables\',which_VAR,'_','table4_rMAEs_h1_final'),title_,subtitle_,1);

table5_rMAEs_h4 = string(sprintfc('%0.2f',  rMAEs_h4 ));
ndx_asterisks = rMAEs_h4<1;
asterisks_rMAEs_h4 = string([ndx_asterisks.*GWtests_all_models_MAE_h4]);
asterisks_rMAEs_h4=strrep(asterisks_rMAEs_h4,'0','');
asterisks_rMAEs_h4=strrep(asterisks_rMAEs_h4,'1','*');
asterisks_rMAEs_h4=strrep(asterisks_rMAEs_h4,'2','**');
asterisks_rMAEs_h4=strrep(asterisks_rMAEs_h4,'3','***');
table5_rMAEs_h4 = cellstr(append(table5_rMAEs_h4,asterisks_rMAEs_h4));
table5_rMAEs_h4_final = array2table([cellstr(all_models_pretty) table5_rMAEs_h4]);
table5_rMAEs_h4_final.Properties.VariableNames = ["Model" pretty_names_short];
table5_rMAEs_h4_final(ndxBENCH,:)=[];
title_=append(title_VAR,'MAE Forecast performance vs. SV-VAR, h = 4');
subtitle_=['Note: Comparison with the "SV-VAR" (Minnesota Stochastic Volatility, baseline ' ...
    'in denominator) against different specifications. Values below 1 ' ...
    'indicate improvement over baseline. Evaluation windows with forecast ' ...
    'origins from 2008Q1 through 2018:Q4. Real time vintages, from the ALFRED database. Significance is assessed by the Giacomini-White test using ' ...
    'Newey-West standard errors. Legend: SV: Stochastic Volatility, SS: Steady State, SR: Shadow Rate'];

table2latex(table5_rMAEs_h4_final,append('results\',which_VAR,'\tables\',which_VAR,'_','table5_rMAEs_h4_final'),title_,subtitle_,1);

table6_CRPS_h1 = string(sprintfc('%0.2f', rCRPS_h1 ));
ndx_asterisks = rCRPS_h1<1;
asterisks_CRPS_h1 = string([ndx_asterisks.*GWtests_all_models_CRPS_h1]);
asterisks_CRPS_h1=strrep(asterisks_CRPS_h1,'0','');
asterisks_CRPS_h1=strrep(asterisks_CRPS_h1,'1','*');
asterisks_CRPS_h1=strrep(asterisks_CRPS_h1,'2','**');
asterisks_CRPS_h1=strrep(asterisks_CRPS_h1,'3','***');
table6_CRPS_h1 = cellstr(append(table6_CRPS_h1,asterisks_CRPS_h1));
table6_CRPS_h1_final = array2table([cellstr(all_models_pretty) table6_CRPS_h1]);
table6_CRPS_h1_final.Properties.VariableNames = ["Model" pretty_names_short];
table6_CRPS_h1_final(ndxBENCH,:)=[];
title_=append(title_VAR,'CRPS Forecast performance vs. SV-VAR, h = 1');
subtitle_=['Note: Comparison with the "SV-VAR" (Minnesota Stochastic Volatility, baseline ' ...
    'in denominator) against different specifications. Values below 1 ' ...
    'indicate improvement over baseline. Evaluation windows with forecast ' ...
    'origins from 2008Q4 through 2019:Q3. Real time vintages, from the ALFRED database. Significance is assessed by the Giacomini-White test using ' ...
    'Newey-West standard errors. Legend: SV: Stochastic Volatility, SS: Steady State, SR: Shadow Rate'];

table2latex(table6_CRPS_h1_final,append('results\',which_VAR,'\tables\',which_VAR,'_','table6_CRPS_h1_final'),title_,subtitle_,1);

table7_CRPS_h4 = string(sprintfc('%0.2f', rCRPS_h4 ));
ndx_asterisks = rCRPS_h4<1;
asterisks_CRPS_h4 = string(ndx_asterisks.*GWtests_all_models_CRPS_h4);
asterisks_CRPS_h4=strrep(asterisks_CRPS_h4,'0','');
asterisks_CRPS_h4=strrep(asterisks_CRPS_h4,'1','*');
asterisks_CRPS_h4=strrep(asterisks_CRPS_h4,'2','**');
asterisks_CRPS_h4=strrep(asterisks_CRPS_h4,'3','***');
table7_CRPS_h4 = cellstr(append(table7_CRPS_h4,asterisks_CRPS_h4));
table7_CRPS_h4_final = array2table([cellstr(all_models_pretty) table7_CRPS_h4]);
table7_CRPS_h4_final.Properties.VariableNames = ["Model" pretty_names_short];
table7_CRPS_h4_final(ndxBENCH,:)=[];

title_=append(title_VAR,'CRPS Forecast performance vs. SV-VAR, h = 4');
subtitle_=['Note: Comparison with the "SV-VAR" (Minnesota Stochastic Volatility, baseline ' ...
    'in denominator) against different specifications. Values below 1 ' ...
    'indicate improvement over baseline. Evaluation windows with forecast ' ...
    'origins from 2008Q1 through 2018:Q4. Real time vintages, from the ALFRED database. Significance is assessed by the Giacomini-White test using ' ...
    'Newey-West standard errors. Legend: SV: Stochastic Volatility, SS: Steady State, SR: Shadow Rate'];

table2latex(table7_CRPS_h4_final,append('results\',which_VAR,'\tables\',which_VAR,'_','table7_CRPS_h4_final'),title_,subtitle_,1);
