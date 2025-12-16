function  [table_results] = create_tables_final(naming,irf_forecast_horizon,no_of_models,draw_GIRFs, ndxBENCH, ...
    all_models_pretty, pretty_names, var_mnemonic_i, tcode, minnesotaPriorMean, which_VAR, results, ...
    tables_dir, GIRFs_dir, shadow_rates_graphs_dir, volatilities_dir, forecast_fancharts, CHOOSE_VAR,...
    IRFPlus_tails_all_models_samples, IRFMinus_tails_all_models_samples, vintages_to_run, mnemonics_for_girfs, ...
    all_models, Yraw_table_last_vintage, shadowrateTails_all_models_samples, FPActuals_irf, ndxMODEL, fcstYdraws_all_models_samples, ...
    post_h_all_models_samples, no_of_samples, spf_dataset_SSP, SSP_id, MU_all_models_samples,forecast_horizon_set,WUXIAshadow2, ...
    produce_tables, filter_variables)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%           TABLES             %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pretty_names_short = char(var_mnemonic_i);
pretty_names_short = string((pretty_names_short(:,1:3,:)));
make_strings_for_graphs

% all_models_pretty = strrep(all_models_pretty,"SV","Stochastic Volatility");
% all_models_pretty = strrep(all_models_pretty,"SR","Shadow Rate");
% all_models_pretty = strrep(all_models_pretty,"SSP","Steady State");
% all_models_pretty = strrep(all_models_pretty,"SStochastic VolatilityS","SSVS");
% all_models_pretty = strrep(all_models_pretty,"SSVS","Stochastic Search Variable Selection");
all_latex_tables=struct();
if produce_tables == 1
    table_list_of_variables = table(string(pretty_names'), string(var_mnemonic_i'), string(pretty_names_short'), string(string(tcode')), string(string(minnesotaPriorMean)));
    table_list_of_variables = renamevars(table_list_of_variables,["Var1","Var2","Var3","Var4","Var5"],["Variable","ALFRED code", "Table mnemonic","Transformation code","Minnesota Prior"]);
    title_ = 'List of Variables';
    subtitle_ = ['Note: Data gathered from the ALFRED-MD. Includes vintages from 2008Q1 to 2019Q4. The transformations legend is the following: 6: ' ...
        'log second differences, 5: log differences, 4: log transformation, 3: second differences, 2: first differences, 1: no transformation. ' ...
        'The Table mnemonic is the shortened name of the variable as used in the results table, while the Minnesota prior column denotes whether ' ...
        'shrinkage of the variable coefficient is towards zero or towards 1 (random walk).'];
    
    filename = char([tables_dir+'table_list_of_variables2.tex']);
    
    % Open the file for writing
    fileID = fopen(filename, 'w');
    
    % Write the LaTeX table header
    fprintf(fileID, '\\begin{table}[h!]\n');
    fprintf(fileID, '\\centering\n');
    fprintf(fileID, '\\begin{tabular}{|l|l|l|c|c|}\n');
    fprintf(fileID, '\\hline\n');
    fprintf(fileID, 'Variable & ALFRED code & Table mnemonic & Transformation code & Minnesota Prior \\\\ \\hline\n');
    
    % Write the table data
    for i = 1:height(table_list_of_variables)
        fprintf(fileID, '%s & %s & %s & %.2f & %.2f \\\\ \\hline\n', ...
            table_list_of_variables.Variable{i}, ...
            table_list_of_variables.('ALFRED code'){i}, ...
            table_list_of_variables.('Table mnemonic'){i}, ...
            table_list_of_variables.('Transformation code')(i), ...
            table_list_of_variables.('Minnesota Prior')(i));
    end
    
    % Write the LaTeX table footer
    fprintf(fileID, '\\end{tabular}\n');
    fprintf(fileID, '\\caption{List of Variables}\n');
    fprintf(fileID, '\\label{tab:variables}\n');
    fprintf(fileID, '\\end{table}\n');
    
    % Close the file
    fclose(fileID);
    
    format bank
    % ndxMAINVARS = find(ismember(var_mnemonic, {'GDP','UNRATE','CPIAUCSL','FEDFUNDS'}));
    table_results=struct();
    subtitle_=[append('Note: Comparison with the ', all_models_pretty(ndxBENCH), [' model (baseline in denominator) against different specifications. ' ...
                'Values below 1 indicate improvement over baseline. Evaluation windows with forecast origins from 2007Q1 through 2019:Q3. ' ...
                'Real time vintages, from the ALFRED database. Significance is assessed by the Giacomini-White test using Newey-West standard errors. ' ...
                'Legend: SV: Stochastic Volatility, SSP: Steady State, SR: Shadow Rate. Highlighted with cyan color are the shadow rate model specifications. ' ...
                'Stars represent the GW test for the rejection of the hypothesis H0: "no significant difference between the forecast of model j vs benchmark", ' ...
                'for the 10\% (*), 5\% (**) and 1\% (***) significance level. Shades from red (benchmark better) to blue (alternative better) show relative performance.'])];
        
    % forecast_horizon_set_graphs = forecast_horizon_set(2:4);
    % fieldNames  = {'h'+string(forecast_horizon_set(1)), 'h'+string(forecast_horizon_set(2)), 'h'+string(forecast_horizon_set(3))};
    fieldNames = cell(1, length(forecast_horizon_set)); % Pre-allocate for efficiency
    for i = 1:length(forecast_horizon_set)
        fieldNames{i} = ['h', num2str(forecast_horizon_set(i))];
    end
    
    subsamples  = {'','_09q1_15q4','_16q1_19q4'};
    scores      = {'RMSE','MAE', 'CRPS'};
    for score_i = scores 
        for subsample_i = subsamples
            for horizon_i = 1:length(forecast_horizon_set)
               
                % Debugging
                % score_i=scores(1);
                % subsample_i=subsamples(1);
                % horizon_i = 1;
                
                score_string = append('r',score_i,'s');
                if string(score_i)=="CRPS"
                    score_string=strrep(score_string,'s','');
                end
                % score table
                table_rSCOREs = string(sprintfc('%0.2f', results.(string(fieldNames(horizon_i))).(string(append(score_string,subsample_i))) ));
                % ndx_asterisks = results.(string(fieldNames(horizon_i))).(string(append(score_string,subsample_i)))<1;
                ndx_asterisks = 1;
                asterisks_rSCOREs = string([ndx_asterisks.*results.(string(fieldNames(horizon_i))).(string(append('GWtests_all_models_',score_i,subsample_i)))]);
                asterisks_rSCOREs = strrep(asterisks_rSCOREs,'0','');
                asterisks_rSCOREs = strrep(asterisks_rSCOREs,'1','*');
                asterisks_rSCOREs = strrep(asterisks_rSCOREs,'2','**');
                asterisks_rSCOREs = strrep(asterisks_rSCOREs,'3','***');
                table_rSCOREs = cellstr(append(table_rSCOREs,asterisks_rSCOREs));
                table_rSCOREs_final = array2table([cellstr(all_models_pretty) table_rSCOREs]);
                table_rSCOREs_final.Properties.VariableNames = ["Model" pretty_names_short];
                table_rSCOREs_final(ndxBENCH,:)=[];
                
                table_rSCOREs_final2_hor{horizon_i} = table_rSCOREs_final(:,[1 filter_variables+1]);
                table_rSCOREs_final2_hor{horizon_i}.Properties.VariableNames = append(string(table_rSCOREs_final2_hor{horizon_i}.Properties.VariableNames), string(fieldNames(horizon_i)));
                table_rSCOREs_final2_hor{horizon_i}.Properties.VariableNames(1)="Model";
            end
            counterr=1;
            finalTable=[];
            for var_ii=2:length(table_rSCOREs_final2_hor{1}.Properties.VariableNames)
                % Combine the tables using outer join
                T_combined = outerjoin(table_rSCOREs_final2_hor{1}(:,[1 var_ii]), table_rSCOREs_final2_hor{2}(:,[1 var_ii]), 'MergeKeys', true);
                if length(fieldNames)>2
                    for kj = 3:length(fieldNames)
                        T_combined = outerjoin(T_combined, table_rSCOREs_final2_hor{kj}(:,[1 var_ii]), 'MergeKeys', true);
                    end
                end
                % T_combined = outerjoin(T_combined, table_rSCOREs_final2_hor{4}(:,[1 var_ii]), 'MergeKeys', true);
                % 
                % Create the final table with 'model' and 'gdp' as subcolumns
                finalTables{counterr} = [T_combined.Model, T_combined(:,2:size(T_combined,2))];
                % finalTables{counterr} = table(T_combined.Model, ...
                %     [table2array(T_combined(:,2:4))], ...
                %     'VariableNames', {'Model', char(T_combined.Properties.VariableNames(2))});
                % 
                if counterr ~= 1
                    finalTable=outerjoin(finalTable,finalTables{counterr}, 'MergeKeys', true);
                else
                    finalTable=finalTables{counterr};
                end
    
                counterr=counterr+1;
                
            end
            % finalTable.Properties.VariableNames=strrep(finalTable.Properties.VariableNames,char(fieldNames{1}),'');
            finalTable.Properties.VariableNames(1)={'Model'};
            [~, idx] = ismember(table_rSCOREs_final2_hor{horizon_i}.Model, finalTable.Model);
            finalTable = finalTable(idx, :);
    
            title_=append(score_i, ' forecast performance vs. , ', all_models_pretty(ndxBENCH), '-VAR, ', char(strrep(strrep(subsample_i,'_',' 20'),'q',' Q')));
            
            latex_file_i = table2latex2color(finalTable,convertStringsToChars(append(char(tables_dir), ...
                'bnch',all_models_pretty(ndxBENCH),string(subsample_i),'_table_r',score_i,'s')),score_i,title_,subtitle_,1,all_models_pretty, fieldNames);

            table_results.(string(append(score_string,subsample_i))) = finalTable;

            all_latex_tables.(string(append(score_string,subsample_i)))=latex_file_i;
        end
        
    end

    table_results.all_latex_tables =all_latex_tables;
end
