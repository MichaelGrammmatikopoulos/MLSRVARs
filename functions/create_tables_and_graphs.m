function  [table_graph_results] = create_tables_and_graphs(naming,irf_forecast_horizon,no_of_models,draw_GIRFs, ndxBENCH,all_models_pretty, pretty_names, var_mnemonic_i, ...
    tcode, minnesotaPriorMean, which_VAR, results, ...
    tables_dir, GIRFs_dir, shadow_rates_graphs_dir, volatilities_dir, forecast_fancharts, CHOOSE_VAR,...
    IRFPlus_tails_all_models_samples, IRFMinus_tails_all_models_samples, vintages_to_run, mnemonics_for_girfs, ...
    all_models, Yraw_table_last_vintage, shadowrateTails_all_models_samples, ...
    FPActuals_irf, ndxMODEL, fcstYdraws_all_models_samples, ...
    post_h_all_models_samples, no_of_samples, spf_dataset_SSP, SSP_id, MU_all_models_samples,forecast_horizon_set,WUXIAshadow2, ...
    produce_SR_graphs, produce_IRFs_graphs, produce_SV_graphs, produce_FAN_graphs, filter_variables)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%           TABLES             %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pretty_names_short = char(var_mnemonic_i);
pretty_names_short = string((pretty_names_short(:,1:3,:)));
make_strings_for_graphs

all_models_pretty = strrep(all_models_pretty,"SV","Stochastic Volatility");
all_models_pretty = strrep(all_models_pretty,"SR","Shadow Rate");
all_models_pretty = strrep(all_models_pretty,"SSP","Steady State");
all_models_pretty = strrep(all_models_pretty,"SStochastic VolatilityS","SSVS");
all_models_pretty = strrep(all_models_pretty,"SSVS","Stochastic Search Variable Selection");

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
table_graph_results=struct();
subtitle_=[append('Note: Comparison with the ', all_models_pretty(ndxBENCH), [' model (baseline in denominator) against different specifications. ' ...
            'Values below 1 indicate improvement over baseline. Evaluation windows with forecast origins from 2007Q1 through 2019:Q3. ' ...
            'Real time vintages, from the ALFRED database. Significance is assessed by the Giacomini-White test using Newey-West standard errors. ' ...
            'Legend: SV: Stochastic Volatility, SSP: Steady State, SR: Shadow Rate. Highlighted with cyan color are the shadow rate model specifications. ' ...
            'Stars represent the GW test for the rejection of the hypothesis H0: "no significant difference between the forecast of model j vs benchmark", for the 10\% (*), 5\% (**) and 1\% (***) significance level.'])];
    
fieldNames  = {'h'+string(forecast_horizon_set(1)), 'h'+string(forecast_horizon_set(2)), 'h'+string(forecast_horizon_set(3))};
subsamples  = {'','_09q1_15q4','_16q1_19q4'};
scores      = {'RMSE','MAE', 'CRPS'};
for score_i = scores 
    for subsample_i = subsamples
        for horizon_i = 1:length(forecast_horizon_set)
           
            score_string = append('r',score_i,'s');
            if string(score_i)=="CRPS"
                score_string=strrep(score_string,'s','');
            end
            % score table
            table_rSCOREs = string(sprintfc('%0.2f', results.(string(fieldNames(horizon_i))).(string(append(score_string,subsample_i))) ));
            ndx_asterisks = results.(string(fieldNames(horizon_i))).(string(append(score_string,subsample_i)))<1;
            asterisks_rSCOREs = string([ndx_asterisks.*results.(string(fieldNames(horizon_i))).(string(append('GWtests_all_models_',score_i,subsample_i)))]);
            asterisks_rSCOREs = strrep(asterisks_rSCOREs,'0','');
            asterisks_rSCOREs = strrep(asterisks_rSCOREs,'1','*');
            asterisks_rSCOREs = strrep(asterisks_rSCOREs,'2','**');
            asterisks_rSCOREs = strrep(asterisks_rSCOREs,'3','***');
            table_rSCOREs = cellstr(append(table_rSCOREs,asterisks_rSCOREs));
            table_rSCOREs_final = array2table([cellstr(all_models_pretty) table_rSCOREs]);
            table_rSCOREs_final.Properties.VariableNames = ["Model" pretty_names_short];
            table_rSCOREs_final(ndxBENCH,:)=[];
            title_=append(score_i, ' forecast performance vs. , ', all_models_pretty(ndxBENCH), '-VAR, ', strrep((string(fieldNames(horizon_i))),'h','h = '), char(strrep(strrep(subsample_i,'_',' 20'),'q',' Q')));
            
            table_rSCOREs_final2 = table_rSCOREs_final(:,[1 filter_variables+1]);

            table2latex(table_rSCOREs_final2,convertStringsToChars(append(char(tables_dir),which_VAR,'_',string(fieldNames(horizon_i)),string(subsample_i),'_table_r',score_i,'s_final')),title_,subtitle_,1,all_models);
  
            table_graph_results.(string(fieldNames(horizon_i))).(string(append(score_string,subsample_i))) = table_rSCOREs_final2;
        end
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%           GRAPHS             %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontsize = 12;
colorPlus     = Colors4Plots(1);
colorMinus    = Colors4Plots(2);
colorBase     = Colors4Plots(8);
colorPlus_vint1     = Colors4Plots(1);
colorMinus_vint1    = Colors4Plots(2);
colorPlus_vint2     = Colors4Plots(2);
colorMinus_vint2    = Colors4Plots(1);



if draw_GIRFs == 1 && produce_IRFs_graphs == 1

    vint_i1 = vintages_to_run(1);
    vint_i2 = vintages_to_run(2);

    if ~exist('title_VAR')
        title_VAR = '';
    end
    %% PLOT GIRF RESULTS
    for model_i = 1:no_of_models

        if CHOOSE_VAR == 1
            figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','GDP','UNRATE','CPIAUCSL', ...
                }));
        elseif CHOOSE_VAR == 2
            figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','GDP','UNRATE','CPIAUCSL', ...
                'GS1', 'GS10'}));
        elseif CHOOSE_VAR == 3
            figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','GDP','UNRATE','CPIAUCSL', ...
                'INDPRO','PAYEMS','PPIACO','PCECC96','HOUST'}));
        elseif CHOOSE_VAR == 4
            figure_varlist = find(ismember(var_mnemonic_i,mnemonics_for_girfs));
        end

        %   figure_varlist=[1:9];
        figure_varnames = [var_mnemonic_i(figure_varlist)];
        figure_pretty_names = [pretty_names(figure_varlist)];

        % mkdir(cd,'results')

        %  for model_i=1:18
        thisfig_plus = figure;
        sgtitle(all_models_pretty(model_i));

        for n = 1:length(figure_varlist)
            if CHOOSE_VAR == 1
                subplot(2,2,n);
            elseif CHOOSE_VAR == 2
                subplot(2,3,n);
            elseif CHOOSE_VAR == 3
                subplot(3,3,n);
            elseif CHOOSE_VAR == 4
                subplot(3,3,n);
            end

            hold on;
            grid on;

            set(gca, 'FontSize', fontsize);
            xaxis = 1:irf_forecast_horizon;
            hplus1  = plot(xaxis, IRFPlus_tails_all_models_samples{model_i}{vint_i1}(figure_varlist(n),:,3), '-.', 'color', colorPlus_vint1, 'linewidth', 2);
            plot(xaxis, squeeze(IRFPlus_tails_all_models_samples{model_i}{vint_i1}(figure_varlist(n),:,:)), '-.', 'color', colorPlus_vint1, 'linewidth', 1);
            hplus2  = plot(xaxis, IRFPlus_tails_all_models_samples{model_i}{vint_i2}(figure_varlist(n),:,3), '-', 'color', colorPlus_vint2, 'linewidth', 2);
            plot(xaxis, squeeze(IRFPlus_tails_all_models_samples{model_i}{vint_i2}(figure_varlist(n),:,:)), '-', 'color', colorPlus_vint2, 'linewidth', 1);
            xlim([1 irf_forecast_horizon]);
            xticks = ([4:4:irf_forecast_horizon]);
            set(gca, 'XTick', xticks);
            yline(0, 'k:');
            %             if n == 3
            %                 lgd = legend([hplus1 hplus2], [naming(vint_i1) naming(vint_i2)],'Location','northeast');
            %             end
            title(figure_pretty_names(n), 'FontWeight', 'normal');

        end
        pos = get(thisfig_plus, 'Position');
        if CHOOSE_VAR == 1
            subplot(2,2,n);
        elseif CHOOSE_VAR == 2
            set(thisfig_plus, 'Position',pos+[-350 -350 350 100]);
        elseif CHOOSE_VAR == 3
            set(thisfig_plus, 'Position',pos+[-350 -350 350 350]);
        elseif CHOOSE_VAR == 4
            set(thisfig_plus, 'Position',pos+[-350 -350 350 350]);
        end

        % Assuming GIRFs_dir and which_VAR are defined
        file_name = append(char(GIRFs_dir), char(append(which_VAR, '_', 'GIRF_plus_', all_models_pretty(model_i), '.png')));
        %         print(append(char(GIRFs_dir), char(append(which_VAR, '_', 'GIRF_plus_', all_models_pretty(model_i)))), '-dpng');

        % Create the legend for the whole figure after all subplots
        % Use 'dummy' plots to create legend entries
        hplus1 = plot(nan, nan, '-.', 'color', colorPlus_vint1, 'linewidth', 2);
        hplus2 = plot(nan, nan, '-', 'color', colorPlus_vint2, 'linewidth', 2);

        % Add the legend outside the subplots
        lgd = legend([hplus1 hplus2], [naming(vint_i1) naming(vint_i2)]);
        lgd.Position = [0.47, 0.02, 0.1, 0.05]; % Set the position of the legend manually
        lgd.Box = 'off'; % Turn off the legend box
        lgd.Orientation = 'horizontal'; % Arrange the legend entries horizontally

        % Adjust the figure size to accommodate the legend
        pos = get(thisfig_plus, 'Position');
        set(thisfig_plus, 'Position', [pos(1), pos(2), pos(3), pos(4)]); % Adjust as needed

        set(gcf, 'Color', 'w')
        export_fig(file_name, '-png', '-a1', '-native')
    end

    for model_i = 1:no_of_models

        if CHOOSE_VAR == 1
            figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','GDP','UNRATE','CPIAUCSL', ...
                }));
        elseif CHOOSE_VAR == 2
            figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','GDP','UNRATE','CPIAUCSL', ...
                'GS1', 'GS10'}));
        elseif CHOOSE_VAR == 3
            figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','GDP','UNRATE','CPIAUCSL', ...
                'INDPRO','PAYEMS','PPIACO','PCECC96','HOUST'}));
        elseif CHOOSE_VAR == 4
            figure_varlist = find(ismember(var_mnemonic_i,mnemonics_for_girfs));
        end

        %   figure_varlist=[1:9];
        figure_varnames = [var_mnemonic_i(figure_varlist)];
        figure_pretty_names = [pretty_names(figure_varlist)];

        % mkdir(cd,'results')

        %  for model_i=1:18
        thisfig_plus = figure;
        sgtitle(all_models_pretty(model_i));

        for n = 1:length(figure_varlist)
            if CHOOSE_VAR == 1
                subplot(2,2,n);
            elseif CHOOSE_VAR == 2
                subplot(2,3,n);
            elseif CHOOSE_VAR == 3
                subplot(3,3,n);
            elseif CHOOSE_VAR == 4
                subplot(3,3,n);
            end

            hold on;
            grid on;

            set(gca, 'FontSize', fontsize);
            xaxis = 1:irf_forecast_horizon;
            hplus1  = plot(xaxis, IRFMinus_tails_all_models_samples{model_i}{vint_i1}(figure_varlist(n),:,3), '-.', 'color', colorPlus_vint1, 'linewidth', 2);
            plot(xaxis, squeeze(IRFMinus_tails_all_models_samples{model_i}{vint_i1}(figure_varlist(n),:,:)), '-.', 'color', colorPlus_vint1, 'linewidth', 1);
            hplus2  = plot(xaxis, IRFMinus_tails_all_models_samples{model_i}{vint_i2}(figure_varlist(n),:,3), '-', 'color', colorPlus_vint2, 'linewidth', 2);
            plot(xaxis, squeeze(IRFMinus_tails_all_models_samples{model_i}{vint_i2}(figure_varlist(n),:,:)), '-', 'color', colorPlus_vint2, 'linewidth', 1);
            xlim([1 irf_forecast_horizon]);
            xticks = ([4:4:irf_forecast_horizon]);
            set(gca, 'XTick', xticks);
            yline(0, 'k:');
            if n == 3
                lgd = legend([hplus1 hplus2], [naming(vint_i1) naming(vint_i2)],'Location','northeast');
            end
            title(figure_pretty_names(n), 'FontWeight', 'normal');

        end
        pos = get(thisfig_plus, 'Position');
        if CHOOSE_VAR == 1
            subplot(2,2,n);
        elseif CHOOSE_VAR == 2
            set(thisfig_plus, 'Position',pos+[-350 -350 350 100]);
        elseif CHOOSE_VAR == 3
            set(thisfig_plus, 'Position',pos+[-350 -350 350 350]);
        elseif CHOOSE_VAR == 4
            set(thisfig_plus, 'Position',pos+[-350 -350 350 350]);
        end
        % Assuming GIRFs_dir and which_VAR are defined
        file_name = append(char(GIRFs_dir), char(append(which_VAR, '_', 'GIRF_minus_', all_models_pretty(model_i))));
        %         print(append(char(GIRFs_dir), char(append(which_VAR, '_', 'GIRF_minus_', all_models_pretty(model_i)))), '-dpng');
        % Create the legend for the whole figure after all subplots
        % Use 'dummy' plots to create legend entries
        hplus1 = plot(nan, nan, '-.', 'color', colorPlus_vint1, 'linewidth', 2);
        hplus2 = plot(nan, nan, '-', 'color', colorPlus_vint2, 'linewidth', 2);

        % Add the legend outside the subplots
        lgd = legend([hplus1 hplus2], [naming(vint_i1) naming(vint_i2)]);
        lgd.Position = [0.47, 0.02, 0.1, 0.05]; % Set the position of the legend manually
        lgd.Box = 'off'; % Turn off the legend box
        lgd.Orientation = 'horizontal'; % Arrange the legend entries horizontally

        % Adjust the figure size to accommodate the legend
        pos = get(thisfig_plus, 'Position');
        set(thisfig_plus, 'Position', [pos(1), pos(2), pos(3), pos(4)]); % Adjust as needed
        set(gcf, 'Color', 'w')
        export_fig(file_name, '-png', '-a1', '-native')
    end
end

    %% ----------------------------------------------------------------------------------------
    % Graph 2: Shadow rate estimates, final and real time
    % -----------------------------------------------------------------------------------------

    % WUXIAshadow = matlab.io.datastore.FileSet("WuXiaShadowRate.xlsx");
    % WUXIAshadow2 = spreadsheetDatastore(WUXIAshadow);
    % WUXIAshadow2.Sheets=1;
    % WUXIAshadow2 = read(WUXIAshadow2);

    if produce_SR_graphs == 1
        corr_tbl = table([all_models_pretty; "Wu-Xia"], zeros(size(all_models_pretty,1)+1,1), zeros(size(all_models_pretty,1)+1,1), ... 
            zeros(size(all_models_pretty,1)+1,1), zeros(size(all_models_pretty,1)+1,1),zeros(size(all_models_pretty,1)+1,1), zeros(size(all_models_pretty,1)+1,1), ... 
            zeros(size(all_models_pretty,1)+1,1), zeros(size(all_models_pretty,1)+1,1));
        corr_tbl_names = ["Model","correlation (real time)", "min (real time)", "mean (real time)","standard deviation (real time)","correlation (full sample)", "min (full sample)", "mean (full sample)","standard deviation (full sample)"];
    corr_tbl.Properties.VariableNames=corr_tbl_names;
        for mm_i = 1:no_of_models
    
            if contains(all_models(mm_i),'srp','ignorecase',true)
                thisfig_shadow = figure;
                sgtitle(append(title_VAR,' - ',"Shadow Rate estimates for model: ", all_models_pretty(mm_i)));
    
                subplot(1,2,1);
    
                hold on
                start_of_graph="1990 Q1";
                last_non_shadow_point="2008 Q4";
                end_of_graph="2019 Q3";
                start_of_graph_num = str2double(extractBetween( start_of_graph , 1 , 4 ));
    
                start_date=find(strcmp(start_of_graph,cell2mat(Yraw_table_last_vintage.observation)));
                last_non_shadow_date=find(strcmp(last_non_shadow_point,cell2mat(Yraw_table_last_vintage.observation)));
                end_date=find(strcmp(end_of_graph,cell2mat(Yraw_table_last_vintage.observation)));
    
                shadow_part_gr = reshape(shadowrateTails_all_models_samples{mm_i}{no_of_samples}(:,1,[1,2,3]),[],3);
                wu_xia_part_gr = table2array(WUXIAshadow2(1:find(strcmp(end_of_graph,cell2mat(WUXIAshadow2.date))),"WuXiaShadow"));
                sh_gr=[Yraw_table_last_vintage.FEDFUNDS(start_date:end_date) [repmat(Yraw_table_last_vintage.FEDFUNDS(start_date:last_non_shadow_date),1,3); shadow_part_gr]];
    
                dn=datenum(start_of_graph_num,1+[0:3:12*(2020-start_of_graph_num)].',1);  % make a sample time vector
                dn=dn(1:end-2);
    
                actuals=plot(dn,sh_gr(:,1),'-', 'color', colorPlus_vint1, 'linewidth', 2);
                sh_gr_bands=plot(dn,sh_gr(:,3),'-', 'color', colorPlus_vint2, 'linewidth', 2);
                plot(dn, squeeze(sh_gr(:,[2,4])), '-.', 'color', colorPlus_vint2, 'linewidth', 1);
                wuxia=plot(dn, wu_xia_part_gr, '-.', 'color', "black", 'linewidth', 2);
    
                idx_sh = sh_gr(:,1)<=0.25;
                corr_tbl(mm_i,6) = table(corr(sh_gr(idx_sh,3),wu_xia_part_gr(idx_sh)));
                corr_tbl(mm_i,7) = table(min(sh_gr(idx_sh,3)));
                corr_tbl(mm_i,8) = table(mean(sh_gr(idx_sh,3)));
                corr_tbl(mm_i,9) = table(std(sh_gr(idx_sh,3)));
    
                xticks = linspace(dn(1), dn(end), 10);  % Adjust the number 10 to change the number of ticks
                set(gca, 'XTick', xticks);
                datetick('x','YYYYQQ','keepticks');
                %             datetick('x','YYYYQQ');        % format axes as time
                xlim([dn(1) dn(end)]);          % fit axes to range of actual data
                lgd = legend([wuxia actuals sh_gr_bands], ["Wu-Xia Shadow rate", "Actual","Full sample estimates (2019Q3 vintage)"],'Location','northeast');
    
                subplot(1,2,2);
    
                hold on
    
                % real time shadow rates
                rt_shadow_rates =[];
                for vint_vint_i=max(forecast_horizon_set)+1:no_of_samples
                    sh_rate_int=reshape(shadowrateTails_all_models_samples{mm_i}{vint_vint_i}(:,1,[1,2,3]),[],3);
                    if vint_vint_i==max(forecast_horizon_set)+1
                        rt_shadow_rates=sh_rate_int(end,:);
                    else
                        rt_shadow_rates=[rt_shadow_rates;sh_rate_int(end,:)];
                    end
                end
    
                shadow_part_gr = rt_shadow_rates;
                wu_xia_part_gr = table2array(WUXIAshadow2(1:find(strcmp(end_of_graph,cell2mat(WUXIAshadow2.date))),"WuXiaShadow"));
                sh_gr=[Yraw_table_last_vintage.FEDFUNDS(start_date:end_date) [repmat(Yraw_table_last_vintage.FEDFUNDS(start_date:last_non_shadow_date),1,3); shadow_part_gr]];
    
                dn=datenum(start_of_graph_num,1+[0:3:12*(2020-start_of_graph_num)].',1);  % make a sample time vector
                dn=dn(1:end-2);
    
                actuals=plot(dn,sh_gr(:,1),'-', 'color', colorPlus_vint1, 'linewidth', 2);
                sh_gr_bands=plot(dn,sh_gr(:,3),'-', 'color', colorPlus_vint2, 'linewidth', 2);
                plot(dn, squeeze(sh_gr(:,[2,4])), '-.', 'color', colorPlus_vint2, 'linewidth', 1);
                wuxia=plot(dn, wu_xia_part_gr, '-.', 'color', "black", 'linewidth', 2);
    
                xticks = linspace(dn(1), dn(end), 10);  % Adjust the number 10 to change the number of ticks
                set(gca, 'XTick', xticks);
                datetick('x','YYYYQQ','keepticks');
                %             datetick('x','YYYYQQ');        % format axes as time
                xlim([dn(1) dn(end)]);          % fit axes to range of actual data
                lgd = legend([wuxia actuals sh_gr_bands], ["Wu-Xia Shadow rate", "Actual","Real time estimates"],'Location','northeast');
    
                pos = get(thisfig_shadow, 'Position');
                set(thisfig_shadow, 'Position',pos+[-100 -100 750 100]);
                file_name = append(char(shadow_rates_graphs_dir), char(append(which_VAR, '_', 'shadow_rate_', all_models_pretty(mm_i))));
                %             print(append(char(shadow_rates_graphs_dir), char(append(which_VAR, '_', 'shadow_rate_', all_models_pretty(mm_i)))), '-dpng');
                set(gcf, 'Color', 'w')
                export_fig(file_name, '-png', '-a1', '-native')
    
                idx_sh = sh_gr(:,1)<=0.25;
                corr_tbl(mm_i,2) = table(corr(sh_gr(idx_sh,3),wu_xia_part_gr(idx_sh)));
                corr_tbl(mm_i,3) = table(min(sh_gr(idx_sh,3)));
                corr_tbl(mm_i,4) = table(mean(sh_gr(idx_sh,3)));
                corr_tbl(mm_i,5) = table(std(sh_gr(idx_sh,3)));
    
            end
    
        end
                corr_tbl(mm_i+1,2) = table(1); 
                corr_tbl(mm_i+1,3) = table(min(wu_xia_part_gr(idx_sh))); 
                corr_tbl(mm_i+1,4) = table(mean(wu_xia_part_gr(idx_sh)));
                corr_tbl(mm_i+1,5) = table(std(wu_xia_part_gr(idx_sh)));
                corr_tbl(mm_i+1,6) = table(1); 
                corr_tbl(mm_i+1,7) = table(min(wu_xia_part_gr(idx_sh))); 
                corr_tbl(mm_i+1,8) = table(mean(wu_xia_part_gr(idx_sh)));
                corr_tbl(mm_i+1,9) = table(std(wu_xia_part_gr(idx_sh)));

    end

    %% ----------------------------------------------------------------------------------------
    % Graph 3: Predictive Distributions
    % -----------------------------------------------------------------------------------------

    if produce_FAN_graphs == 1
        for mmm_i=1:no_of_models
            if CHOOSE_VAR == 1
                figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','GDP','UNRATE','CPIAUCSL'}));
            elseif CHOOSE_VAR == 2
                figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','GDP','UNRATE','CPIAUCSL'}));
            elseif CHOOSE_VAR == 3
                figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','GDP','UNRATE','CPIAUCSL','INDPRO','PAYEMS','PINCOME','PCECC96','HOUST'}));
            elseif CHOOSE_VAR == 4
                %             figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','GDP','UNRATE','CPIAUCSL','INDPRO','PAYEMS','PINCOME','PCECC96','GS10'}));
                figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','TB3MS','TB6MS','GS1','GS3','GS5','GS10','BAA','INFEXP'}));
            end
    
            figure_pretty_names = [pretty_names(figure_varlist)];
    
            this_figure_fanchart=figure;
            set(this_figure_fanchart, 'Position', [100, 100, 1200, 800]); % Adjusting size
    
            FPActuals_irf2=FPActuals_irf(:,ndxMODEL);
            for n = 1:length(figure_varlist)
    
                var_i=figure_varlist(n);
    
                start_of_graph="2008 Q4";
                end_of_graph="2019 Q4";
                start_of_graph_num = str2double(extractBetween( start_of_graph , 1 , 4 ));
    
                hold on
                historical = FPActuals_irf2(:,var_i);
                forecast_rt =[];
                for vint_vint_i=1:no_of_samples
                    irf_int_int=reshape(fcstYdraws_all_models_samples{mmm_i}{vint_vint_i}(var_i,1,:),1,[]);
                    if vint_vint_i==1
                        forecast_rt=irf_int_int;
                    else
                        forecast_rt=[forecast_rt;irf_int_int];
                    end
                end
    
                if CHOOSE_VAR == 1
                    subplot(2,2,n);
                elseif CHOOSE_VAR == 2
                    subplot(2,2,n);
                elseif CHOOSE_VAR == 3
                    subplot(3,3,n);
                elseif CHOOSE_VAR == 4
                    subplot(3,3,n);
                end
                hold on
                sgtitle(append("Predictive distributions vs Actuals (first prints): ", all_models_pretty(mmm_i)));
                [lineh, bandsh] = fanChart(1:size(forecast_rt,1), forecast_rt);
                plot(historical,'LineStyle','--','Marker','o','color','black','MarkerSize',5,'MarkerFaceColor','blue','MarkerEdgeColor','blue','LineWidth',1.2)
                txt = strcat({'Pct'}, cellstr(int2str((10:10:90)')));
    
                title(figure_pretty_names(n), 'FontWeight', 'normal');
    
                if n > length(figure_varlist) - 3 % For last three graphs
                    xticks = (7:8:no_of_samples); % Adjusted for 47 periods
                    xticklabels({'2009 Q4','2011 Q4','2013 Q4','2015 Q4','2017 Q4','2019 Q4'}); % Adjusted for 47 periods
    
                else
                    xticklabels([]); % Removing xlim numbers for other graphs
                end
            end
    
            file_name = append(char(forecast_fancharts), char(append(which_VAR, '_', 'forecast_fancharts_', all_models_pretty(mmm_i))));
            %         print(append(char(forecast_fancharts), char(append(which_VAR, '_', 'forecast_fancharts_', all_models_pretty(mmm_i)))), '-dpng');
            set(gcf, 'Color', 'w')
            export_fig(file_name, '-png', '-a1', '-native')
        end
    end

    %% -----------------------------------------------------------------------------------------
    % Graph 4: Volatilities
    % -----------------------------------------------------------------------------------------
    
    if produce_SV_graphs == 1
        for mmm_i=1:no_of_models
    
            if contains(all_models(mmm_i),'stvol','ignorecase',true)
                if CHOOSE_VAR == 1
                    figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','GDP','UNRATE','CPIAUCSL', ...
                        }));
                elseif CHOOSE_VAR == 2
                    figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','GDP','UNRATE','CPIAUCSL'}));
                elseif CHOOSE_VAR == 3
                    figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','GDP','UNRATE','CPIAUCSL', ...
                        'INDPRO','PAYEMS','PINCOME','PCECC96','HOUST'}));
                elseif CHOOSE_VAR == 4
                    figure_varlist = find(ismember(var_mnemonic_i,{'FEDFUNDS','GDP','UNRATE','CPIAUCSL', ...
                        'INDPRO','PAYEMS','PPIACO','PCECC96','M2SL'}));
                end
    
                figure_pretty_names = [pretty_names(figure_varlist)];
    
                this_figure_volatility_chart=figure;
                for n = 1:length(figure_varlist)
    
                    var_i=figure_varlist(n);
    
                    hold on
                    if CHOOSE_VAR == 1
                        subplot(2,2,n);
                    elseif CHOOSE_VAR == 2
                        subplot(2,2,n);
                    elseif CHOOSE_VAR == 3
                        subplot(3,3,n);
                    elseif CHOOSE_VAR == 4
                        subplot(3,3,n);
                    end
    
                    start_of_graph="1973 Q4";
                    end_of_graph="2019 Q3";
                    start_of_graph_num = str2double(extractBetween( start_of_graph , 1 , 4 ));
                    dn=datenum(start_of_graph_num,1+[0:3:12*(2020-start_of_graph_num)].',1);  % make a sample time vector
                    dn=dn(1:end-2);
    
                    hold on
                    sgtitle(append("Volatility estimates: ", all_models_pretty(mmm_i)));
                    plot(dn,exp(post_h_all_models_samples{mmm_i}{no_of_samples}(:,n)./2))
                    title(figure_pretty_names(n), 'FontWeight', 'normal');
                    datetick('x','YYYYQQ');        % format axes as time
                    xlim([dn(1) dn(end)]);          % fit axes to range of actual data
    
                    % Define the x-axis with 4-6 tick marks with dates
                    ticks = linspace(dn(1), dn(end), 6);  % adjust '6' to the number of ticks you want
                    labels = cellstr(datestr(ticks, 'YYYYQQ'));
                    xticks = (ticks);
                    xticklabels(labels);
                end
                pos = get(this_figure_volatility_chart, 'Position');
                set(this_figure_volatility_chart, 'Position', [pos(1) pos(2) pos(3)+100 pos(4)+100]);
    
                %             print(append(char(volatilities_dir), char(append(which_VAR, '_', 'volatilities_', all_models_pretty(mmm_i)))), '-dpng');
                file_name = append(char(volatilities_dir), char(append(which_VAR, '_', 'volatilities_', all_models_pretty(mmm_i))));
                set(gcf, 'Color', 'w')
                file_name = append(char(volatilities_dir), char(append(which_VAR, '_', 'volatilities_', all_models_pretty(mmm_i))));
                export_fig(file_name, '-png', '-a1', '-native')
            end
        end
    end
    
    % for mmm_i= no_of_models:no_of_models
    %     if contains(all_models(mmm_i),'_ssp','ignorecase',true)
    %         for var_i = 1:21;
    % 
    %             ssp_prior_dt = (table2array(spf_dataset_SSP(:,1+ndxMODEL))'.*SSP_id)';
    % 
    %             thisfig_ssp = figure;
    %             sgtitle(append(title_VAR,' - ',"SSP estimates for model: ", var_mnemonic_i(var_i),all_models_pretty(mmm_i)));
    % 
    %             %         subplot(2,2,1);
    % 
    %             hold on
    % 
    %             % real time SSP
    %             rt_ssp =[];
    % 
    %             %             ssp_int=reshape(MU_all_models_samples(var_i,:,:,model_i),3,no_of_samples)';
    %             ssp_int=cell2mat(cellfun(@(M) M(var_i,:), MU_all_models_samples{mmm_i}, 'UniformOutput', false));
    % 
    %             FPActuals_irf2=FPActuals_irf(:,ndxMODEL);
    %             ssp_part_gr = ssp_int;
    %             historical = FPActuals_irf2(:,var_i);
    % 
    %             start_of_graph="2008 Q2";
    %             end_of_graph="2019 Q4";
    %             start_of_graph_num = str2double(extractBetween( start_of_graph , 1 , 4 ));
    % 
    %             dn=datenum(start_of_graph_num,1+[0:3:12*(2020-start_of_graph_num)].',1);  % make a sample time vector
    %             dn=dn(1:end-2);
    % 
    %             actuals=plot(dn,historical(:,1),'-', 'color', colorPlus_vint1, 'linewidth', 2);
    %             ssp_gr_bands=plot(dn,ssp_part_gr(:,2),'-', 'color', colorPlus_vint2, 'linewidth', 2);
    %             plot(dn, squeeze(ssp_part_gr(:,[1,3])), '-.', 'color', colorPlus_vint2, 'linewidth', 1);
    %             ssp_prior_gr=plot(dn, ssp_prior_dt(:,var_i), '-.', 'color', "black", 'linewidth', 2);
    % 
    %             datetick('x','YYYYQQ');        % format axes as time
    %             xlim([dn(1) dn(end)]);          % fit axes to range of actual data
    %             lgd = legend([ssp_prior_gr ssp_gr_bands actuals ], ["Prior (SPF)","Posterior estimates","Actuals (first prints)"],'Location','northeast');
    % 
    %             %         pos = get(thisfig_shadow, 'Position');
    %             %         set(thisfig_shadow, 'Position',pos+[-100 -100 750 100]);
    % 
    % 
    %             %                 print(append(char(ssp_graphs_dir), char(append(which_VAR, '_', 'volatilities_', all_models_pretty(mmm_i)))), '-dpng');
    % 
    %         end
    %     end
    % end


