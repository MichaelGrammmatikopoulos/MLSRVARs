function  [results] = create_results(yforecast_all_models_samples, ...
    fcstYdraws_all_models_samples, FPActuals_irf, N, no_of_models, ndxMODEL, ...
    no_of_samples, ndxBENCH, forecast_horizon_set, all_forecast_samples_idx,quarters_forecast,window_i)

quarters_forecast = quarters_forecast(max(forecast_horizon_set):end);

idx2009q1 = find(ismember(string(quarters_forecast),{'2009 Q1'}));
idx2015q4 = find(ismember(string(quarters_forecast),{'2015 Q4'}));

idx2016q1 = find(ismember(string(quarters_forecast),{'2016 Q1'}));
idx2019q4 = find(ismember(string(quarters_forecast),{'2019 Q4'}));

idxFE2015q4 = find(ismember(string(quarters_forecast),{'2015 Q4'}));

results = struct();
fieldNames = {'h'+string(forecast_horizon_set(1)), 'h'+string(forecast_horizon_set(2)), 'h'+string(forecast_horizon_set(3)), 'h'+string(forecast_horizon_set(4)), 'h'+string(forecast_horizon_set(5))};
for horizon_i = 1:length(forecast_horizon_set)

    forecast_horizon_i = forecast_horizon_set(horizon_i);
    range_i = all_forecast_samples_idx{horizon_i}(1):all_forecast_samples_idx{horizon_i}(2);

    rRMSEs              = {};
    rRMSEs_09q1_15q4    = {};
    rRMSEs_16q1_19q4    = {};

    rMAEs               = {};
    rMAEs_09q1_15q4     = {};
    rMAEs_16q1_19q4     = {};

    % Collect
    SE_matrices                             = {};
    AE_matrices                             = {};
    CR_matrices                             = {};
    RMSE_scores                           = {};
    RMSE_scores_09q1_15q4                 = {};
    RMSE_scores_16q1_19q4                 = {};
    CRPS_scores                           = {};
    CRPS_scores_09q1_15q4                 = {};
    CRPS_scores_16q1_19q4                 = {};
    MAE_scores                            = {};
    MAE_scores_09q1_15q4                  = {};
    MAE_scores_16q1_19q4                  = {};
    full_forecast_matrices                  = {};

    % Collect ratios
    rCRPS               = {};
    rCRPS_09q1_15q4     = {};
    rCRPS_16q1_19q4     = {};

    for model_i = 1:no_of_models

        full_forecast_matrix = cell2mat(yforecast_all_models_samples{model_i,1});
        full_forecast_matrix = full_forecast_matrix(:,horizon_i);
        full_forecast_matrix = reshape(full_forecast_matrix(:,1,:),N,[])';

        % full_forecast_matrix = full_forecast_matrix(range_i,:);
        full_forecast_matrices{model_i} = full_forecast_matrix;

        % Quantities for the GW test
        SE_matrix  = (FPActuals_irf(range_i,ndxMODEL)-full_forecast_matrix(range_i,:)).^2;
        AE_matrix   = abs(FPActuals_irf(range_i,ndxMODEL)-full_forecast_matrix(range_i,:));

        % % Quantities for the ratios
        RMSE_scores{model_i}           = sqrt(mean(SE_matrix(idx2009q1:idx2015q4,:)));
        RMSE_scores_09q1_15q4{model_i} = sqrt(mean(SE_matrix(idx2009q1:idx2015q4,:)));
        RMSE_scores_16q1_19q4{model_i} = sqrt(mean(SE_matrix(idx2016q1:idx2019q4,:)));
        SE_matrices{model_i} = SE_matrix;

        MAE_scores{model_i}            = mean(AE_matrix(idx2009q1:idx2015q4,:));
        MAE_scores_09q1_15q4{model_i}  = mean(AE_matrix(idx2009q1:idx2015q4,:));
        MAE_scores_16q1_19q4{model_i}  = mean(AE_matrix(idx2016q1:idx2019q4,:));
        AE_matrices{model_i} = AE_matrix;

        % fcstYdraws_all_samples          = cell(no_of_samples,1);
        % fcstYdraws_all_models_samples_i = fcstYdraws_all_models_samples{model_i,1};
        % for sample_i = 1:no_of_samples
        %     fcstYdraws_all_samples{sample_i,1} = fcstYdraws_all_models_samples_i{sample_i,1};
        % end

        %% CRPS
        actuals_i=0;
        for vint_i=1:size(FPActuals_irf,1)
            actuals_i=actuals_i+1;
            for nnn = 1 : N % loop over elements of Y
                CRPSscore(nnn,actuals_i) = crpsDraws(FPActuals_irf(actuals_i,nnn), reshape(fcstYdraws_all_models_samples{model_i}{vint_i}(nnn,horizon_i,:),[],1),[]);
            end
        end
        CR_matrices{model_i} = CRPSscore';
        CR_matrices_09q1_15q4{model_i} = CRPSscore(:,idx2009q1:idx2015q4)';
        CR_matrices_16q1_19q4{model_i} = CRPSscore(:,idx2016q1:idx2019q4)';
                
        CRPS_scores{model_i} = mean(CRPSscore(:,idx2009q1:idx2019q4)');
        CRPS_scores_09q1_15q4{model_i} = mean(CRPSscore(:,idx2009q1:idx2015q4)');
        CRPS_scores_16q1_19q4{model_i} = mean(CRPSscore(:,idx2016q1:idx2019q4)');

    end %for model_i
    
    for model_i = 1:no_of_models
        rRMSEs{model_i}             = round(RMSE_scores{model_i}./RMSE_scores{ndxBENCH},2);
        rRMSEs_09q1_15q4{model_i}   = round(RMSE_scores_09q1_15q4{model_i}./RMSE_scores_09q1_15q4{ndxBENCH},2);
        rRMSEs_16q1_19q4{model_i}   = round(RMSE_scores_16q1_19q4{model_i}./RMSE_scores_16q1_19q4{ndxBENCH},2);

        rMAEs{model_i}              = round(MAE_scores{model_i}./MAE_scores{ndxBENCH},2);
        rMAEs_09q1_15q4{model_i}    = round(MAE_scores_09q1_15q4{model_i}./MAE_scores_09q1_15q4{ndxBENCH},2);
        rMAEs_16q1_19q4{model_i}    = round(MAE_scores_16q1_19q4{model_i}./MAE_scores_16q1_19q4{ndxBENCH},2);

        rCRPS{model_i}              = round(CRPS_scores{model_i}./CRPS_scores{ndxBENCH},2);
        rCRPS_09q1_15q4{model_i}    = round(CRPS_scores_09q1_15q4{model_i}./CRPS_scores_09q1_15q4{ndxBENCH},2);
        rCRPS_16q1_19q4{model_i}    = round(CRPS_scores_16q1_19q4{model_i}./CRPS_scores_16q1_19q4{ndxBENCH},2);

        %% Perform the Giacomini-White test
        if model_i~=ndxBENCH
            for col_i = 1:N
                GWtests_all_models_RMSE(model_i,col_i)           = CPAtest(SE_matrices{model_i}(idx2009q1:idx2019q4,col_i),SE_matrices{ndxBENCH}(idx2009q1:idx2019q4,col_i),forecast_horizon_i,0.05,1);
                GWtests_all_models_RMSE_09q1_15q4(model_i,col_i) = CPAtest(SE_matrices{model_i}(idx2009q1:idx2015q4,col_i),SE_matrices{ndxBENCH}(idx2009q1:idx2015q4,col_i),forecast_horizon_i,0.05,1);
                GWtests_all_models_RMSE_16q1_19q4(model_i,col_i) = CPAtest(SE_matrices{model_i}(idx2016q1:idx2019q4,col_i),SE_matrices{ndxBENCH}(idx2016q1:idx2019q4,col_i),forecast_horizon_i,0.05,1);

                GWtests_all_models_MAE(model_i,col_i)            = CPAtest(AE_matrices{model_i}(idx2009q1:idx2019q4,col_i),AE_matrices{ndxBENCH}(idx2009q1:idx2019q4,col_i),forecast_horizon_i,0.05,1);
                GWtests_all_models_MAE_09q1_15q4(model_i,col_i)  = CPAtest(AE_matrices{model_i}(idx2009q1:idx2015q4,col_i),AE_matrices{ndxBENCH}(idx2009q1:idx2015q4,col_i),forecast_horizon_i,0.05,1);
                GWtests_all_models_MAE_16q1_19q4(model_i,col_i)  = CPAtest(AE_matrices{model_i}(idx2016q1:idx2019q4,col_i),AE_matrices{ndxBENCH}(idx2016q1:idx2019q4,col_i),forecast_horizon_i,0.05,1);

                GWtests_all_models_CRPS(model_i,col_i)           = CPAtest(CR_matrices{model_i}(idx2009q1:idx2019q4,col_i), CR_matrices{ndxBENCH}(idx2009q1:idx2019q4,col_i),forecast_horizon_i,0.05,1);
                GWtests_all_models_CRPS_09q1_15q4(model_i,col_i) = CPAtest(CR_matrices{model_i}(idx2009q1:idx2015q4,col_i), CR_matrices{ndxBENCH}(idx2009q1:idx2015q4,col_i),forecast_horizon_i,0.05,1);
                GWtests_all_models_CRPS_16q1_19q4(model_i,col_i) = CPAtest(CR_matrices{model_i}(idx2016q1:idx2019q4,col_i), CR_matrices{ndxBENCH}(idx2016q1:idx2019q4,col_i),forecast_horizon_i,0.05,1);

            end
        end
    end %for model_i

    results.(fieldNames{horizon_i}).rRMSEs                            = cell2mat(rRMSEs');
    results.(fieldNames{horizon_i}).rRMSEs_09q1_15q4                  = cell2mat(rRMSEs_09q1_15q4');
    results.(fieldNames{horizon_i}).rRMSEs_16q1_19q4                  = cell2mat(rRMSEs_16q1_19q4');
     results.(fieldNames{horizon_i}).RMSE_scores                      = RMSE_scores;
    results.(fieldNames{horizon_i}).SE_matrices                       = SE_matrices; 

    results.(fieldNames{horizon_i}).rMAEs                             = cell2mat(rMAEs');
    results.(fieldNames{horizon_i}).rMAEs_09q1_15q4                   = cell2mat(rMAEs_09q1_15q4');
    results.(fieldNames{horizon_i}).rMAEs_16q1_19q4                   = cell2mat(rMAEs_16q1_19q4');
    results.(fieldNames{horizon_i}).MAE_scores                        = MAE_scores;
    results.(fieldNames{horizon_i}).AE_matrices                       = AE_matrices; 

    % results.(fieldNames{horizon_i}).cRMSEs                            = cRMSEs;
    % results.(fieldNames{horizon_i}).cMAEs                             = cMAEs;
    % results.(fieldNames{horizon_i}).cCRPS                             = cCRPS;

    results.(fieldNames{horizon_i}).rCRPS                             = cell2mat(rCRPS');
    results.(fieldNames{horizon_i}).rCRPS_09q1_15q4                   = cell2mat(rCRPS_09q1_15q4');
    results.(fieldNames{horizon_i}).rCRPS_16q1_19q4                   = cell2mat(rCRPS_16q1_19q4');
    results.(fieldNames{horizon_i}).CRPS_scores                       = CRPS_scores;
    results.(fieldNames{horizon_i}).CR_matrices                       = CR_matrices;

    results.(fieldNames{horizon_i}).GWtests_all_models_RMSE           = GWtests_all_models_RMSE;
    results.(fieldNames{horizon_i}).GWtests_all_models_RMSE_09q1_15q4 = GWtests_all_models_RMSE_09q1_15q4;
    results.(fieldNames{horizon_i}).GWtests_all_models_RMSE_16q1_19q4 = GWtests_all_models_RMSE_16q1_19q4;

    results.(fieldNames{horizon_i}).GWtests_all_models_MAE             = GWtests_all_models_MAE;
    results.(fieldNames{horizon_i}).GWtests_all_models_MAE_09q1_15q4   = GWtests_all_models_MAE_09q1_15q4;
    results.(fieldNames{horizon_i}).GWtests_all_models_MAE_16q1_19q4   = GWtests_all_models_MAE_16q1_19q4;

    results.(fieldNames{horizon_i}).GWtests_all_models_CRPS            = GWtests_all_models_CRPS;
    results.(fieldNames{horizon_i}).GWtests_all_models_CRPS_09q1_15q4  = GWtests_all_models_CRPS_09q1_15q4;
    results.(fieldNames{horizon_i}).GWtests_all_models_CRPS_16q1_19q4  = GWtests_all_models_CRPS_16q1_19q4;

end

% function  [results] = create_results(yforecast_all_models_samples, ...
%     fcstYdraws_all_models_samples, FPActuals, N, no_of_models, ndxMODEL, ...
%     no_of_samples, ndxBENCH, forecast_horizon_set, all_forecast_samples_idx,quarter_names_FEsample,RM_windows)
% 
% % idx2009q1 = find(ismember(string(quarter_names_final),{'2009 Q4'}));
% idx2015q4 = find(ismember(string(quarter_names_FEsample),{'2015 Q4'}));
% idx2019q4 = find(ismember(string(quarter_names_FEsample),{'2019 Q4'}));
% 
% idxFE2015q4 = find(ismember(string(quarter_names_FEsample),{'2015 Q4'}));
% 
% results = struct();
% fieldNames = {'h'+string(forecast_horizon_set(1)), 'h'+string(forecast_horizon_set(2)), 'h'+string(forecast_horizon_set(3)), 'h'+string(forecast_horizon_set(4))};
% for horizon_i = 1:length(forecast_horizon_set)
% 
%     forecast_horizon_i = forecast_horizon_set(horizon_i);
%     range_i = all_forecast_samples_idx{horizon_i}(1):all_forecast_samples_idx{horizon_i}(2);
% 
%     rRMSEs              = {};
%     rRMSEs_09q1_15q4    = {};
%     rRMSEs_16q1_19q4    = {};
% 
%     rMAEs               = {};
%     rMAEs_09q1_15q4     = {};
%     rMAEs_16q1_19q4     = {};
% 
%     % Collect
%     RMSE_scores                           = {};
%     RMSE_scores_09q1_15q4                 = {};
%     RMSE_scores_16q1_19q4                 = {};
%     MAE_scores                            = {};
%     MAE_scores_09q1_15q4                  = {};
%     MAE_scores_16q1_19q4                  = {};
%     CRPS_scores                           = {};
%     CRPS_scores_09q1_15q4                 = {};
%     CRPS_scores_16q1_19q4                 = {};
%     full_forecast_matrices                  = {};
% 
%     % Collect cumulatives: cRMSEs, cMAEs
%     cRMSEs              = NaN(size(FPActuals(:,ndxMODEL),1),size(FPActuals(:,ndxMODEL),2),no_of_models);
%     cMAEs               = NaN(size(FPActuals(:,ndxMODEL),1),size(FPActuals(:,ndxMODEL),2),no_of_models);
%     cCRPS               = NaN(size(FPActuals(:,ndxMODEL),1),size(FPActuals(:,ndxMODEL),2),no_of_models);
% 
%     % Collect ratios
%     rCRPS               = {};
%     rCRPS_09q1_15q4     = {};
%     rCRPS_16q1_19q4     = {};
% 
%     for model_i = 1:no_of_models
% 
%         full_forecast_matrix = cell2mat(yforecast_all_models_samples{model_i,1});
%         full_forecast_matrix = full_forecast_matrix(:,horizon_i);
%         full_forecast_matrix = reshape(full_forecast_matrix(:,1,:),N,[])';
% 
%         full_forecast_matrix = full_forecast_matrix(range_i,:);
%         full_forecast_matrices{model_i} = full_forecast_matrix;
% 
%         % Quantities for the GW test
%         SE_matrix  = (FPActuals(:,ndxMODEL)-full_forecast_matrix).^2;
%         AE_matrix   = abs(FPActuals(:,ndxMODEL)-full_forecast_matrix);
% 
%         % Quantities for the ratios
%         RMSE_matrix            = sqrt(mean((FPActuals(:,ndxMODEL)-full_forecast_matrix).^2));
%         RMSE_matrix_09q1_15q4  = sqrt(mean((FPActuals(1:idx2015q4,ndxMODEL)-full_forecast_matrix(1:idx2015q4,:)).^2));
%         RMSE_matrix_16q1_19q4  = sqrt(mean((FPActuals(idx2015q4+1:end,ndxMODEL)-full_forecast_matrix(idx2015q4+1:end,:)).^2));
% 
%         MAE_matrix             = mean(abs(FPActuals(:,ndxMODEL)-full_forecast_matrix));
%         MAE_matrix_09q1_15q4   = mean(abs(FPActuals(1:idx2015q4,ndxMODEL)-full_forecast_matrix(1:idx2015q4,:)));
%         MAE_matrix_16q1_19q4   = mean(abs(FPActuals(idx2015q4+1:end,ndxMODEL)-full_forecast_matrix(idx2015q4+1:end,:)));
% 
%         % save all quantities
%         RMSE_scores{model_i}           = RMSE_matrix;
%         RMSE_scores_09q1_15q4{model_i} = RMSE_matrix_09q1_15q4;
%         RMSE_scores_16q1_19q4{model_i} = RMSE_matrix_16q1_19q4;
% 
%         MAE_scores{model_i}            = MAE_matrix;
%         MAE_scores_09q1_15q4{model_i}  = MAE_matrix_09q1_15q4;
%         MAE_scores_16q1_19q4{model_i}  = MAE_matrix_16q1_19q4;
% 
%         SE_matrices{model_i} = SE_matrix;
%         AE_matrices{model_i} = AE_matrix;
% 
%         fcstYdraws_all_samples          = cell(no_of_samples,1);
%         fcstYdraws_all_models_samples_i = fcstYdraws_all_models_samples{model_i,1};
%         for sample_i = 1:no_of_samples
%             fcstYdraws_all_samples{sample_i,1} = fcstYdraws_all_models_samples_i{sample_i,1};
%         end
% 
%         %% CRPS
%         actuals_i=0;
%         for vint_i=range_i
%             actuals_i=actuals_i+1;
%             for nnn = 1 : N % loop over elements of Y
%                 CRPSscore(nnn,actuals_i) = crpsDraws(FPActuals(actuals_i,nnn), reshape(fcstYdraws_all_models_samples{model_i}{vint_i}(nnn,horizon_i,:),[],1),[]);
%             end
%         end
% 
%         CR_matrices{model_i} = CRPSscore';
% 
%         CRPS_scores{model_i} = mean(CRPSscore');
%         CRPS_scores_09q1_15q4{model_i} = mean(CRPSscore(:,1:idxFE2015q4)');
%         CRPS_scores_16q1_19q4{model_i} = mean(CRPSscore(:,idxFE2015q4+1:end)');
% 
%     end
% 
%     for model_i = 1:no_of_models
% 
%         rRMSEs{model_i}             = round(RMSE_scores{model_i}./RMSE_scores{ndxBENCH},2);
%         rRMSEs_09q1_15q4{model_i}   = round(RMSE_scores_09q1_15q4{model_i}./RMSE_scores_09q1_15q4{ndxBENCH},2);
%         rRMSEs_16q1_19q4{model_i}   = round(RMSE_scores_16q1_19q4{model_i}./RMSE_scores_16q1_19q4{ndxBENCH},2);
% 
%         rMAEs{model_i}              = round(MAE_scores{model_i}./MAE_scores{ndxBENCH},2);
%         rMAEs_09q1_15q4{model_i}    = round(MAE_scores_09q1_15q4{model_i}./MAE_scores_09q1_15q4{ndxBENCH},2);
%         rMAEs_16q1_19q4{model_i}    = round(MAE_scores_16q1_19q4{model_i}./MAE_scores_16q1_19q4{ndxBENCH},2);
% 
%         rCRPS{model_i}              = round(CRPS_scores{model_i}./CRPS_scores{ndxBENCH},2);
%         rCRPS_09q1_15q4{model_i}    = round(CRPS_scores_09q1_15q4{model_i}./CRPS_scores_09q1_15q4{ndxBENCH},2);
%         rCRPS_16q1_19q4{model_i}    = round(CRPS_scores_16q1_19q4{model_i}./CRPS_scores_16q1_19q4{ndxBENCH},2);
% 
%         cRMSEs{model_i}   = cumsum(SE_matrices{model_i});
%         cMAEs{model_i}    = cumsum(AE_matrices{model_i});
%         cCRPS{model_i}    = cumsum(CR_matrices{model_i});
% 
%         %% Perform the Giacomini-White test
%         if model_i~=ndxBENCH
%             for col_i = 1:N
%                 GWtests_all_models_RMSE(model_i,col_i)           = CPAtest(SE_matrices(:,col_i,model_i),SE_matrices(:,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
%                 GWtests_all_models_RMSE_09q1_15q4(model_i,col_i) = CPAtest(SE_matrices(1:idx2015q4,col_i,model_i),SE_matrices(1:idx2015q4,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
%                 GWtests_all_models_RMSE_16q1_19q4(model_i,col_i) = CPAtest(SE_matrices(idx2015q4+1:end,col_i,model_i),SE_matrices(idx2015q4+1:end,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
% 
%                 GWtests_all_models_MAE(model_i,col_i)            = CPAtest(AE_matrices(:,col_i,model_i),AE_matrices(:,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
%                 GWtests_all_models_MAE_09q1_15q4(model_i,col_i)  = CPAtest(AE_matrices(1:idx2015q4,col_i,model_i),AE_matrices(1:idx2015q4,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
%                 GWtests_all_models_MAE_16q1_19q4(model_i,col_i)  = CPAtest(AE_matrices(idx2015q4+1:end,col_i,model_i),AE_matrices(idx2015q4+1:end,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
% 
%                 GWtests_all_models_CRPS(model_i,col_i)           = CPAtest(CR_matrices(:,col_i,model_i), CR_matrices(:,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
%                 GWtests_all_models_CRPS_09q1_15q4(model_i,col_i) = CPAtest(CR_matrices(1:idx2015q4,col_i,model_i), CR_matrices(1:idx2015q4,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
%                 GWtests_all_models_CRPS_16q1_19q4(model_i,col_i) = CPAtest(CR_matrices(idx2015q4+1:end,col_i,model_i), CR_matrices(idx2015q4+1:end,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
% 
%             end
%         end
%     end
% 
%     results.(fieldNames{horizon_i}).rRMSEs_09q1_15q4                  = reshape(rRMSEs_09q1_15q4,N,no_of_models)';
%     results.(fieldNames{horizon_i}).rRMSEs_16q1_19q4                  = reshape(rRMSEs_16q1_19q4,N,no_of_models)';
%     results.(fieldNames{horizon_i}).rRMSEs                            = reshape(rRMSEs,N,no_of_models)';
%     results.(fieldNames{horizon_i}).RMSE_scores                     = RMSE_scores;
%     results.(fieldNames{horizon_i}).SE_matrices                       = SE_matrices; 
% 
%     results.(fieldNames{horizon_i}).rMAEs_09q1_15q4                   = reshape(rMAEs_09q1_15q4,N,no_of_models)';
%     results.(fieldNames{horizon_i}).rMAEs_16q1_19q4                   = reshape(rMAEs_16q1_19q4,N,no_of_models)';
%     results.(fieldNames{horizon_i}).rMAEs                             = reshape(rMAEs,N,no_of_models)';
%     results.(fieldNames{horizon_i}).MAE_scores                      = MAE_scores;
%     results.(fieldNames{horizon_i}).AE_matrices                       = AE_matrices;
% 
%     results.(fieldNames{horizon_i}).rCRPS                             = reshape(rCRPS,N,no_of_models)';
%     results.(fieldNames{horizon_i}).rCRPS_09q1_15q4                   = reshape(rCRPS_09q1_15q4,N,no_of_models)';
%     results.(fieldNames{horizon_i}).rCRPS_16q1_19q4                   = reshape(rCRPS_16q1_19q4,N,no_of_models)';
%     results.(fieldNames{horizon_i}).CRPS_scores                     = CRPS_scores;
%     results.(fieldNames{horizon_i}).CR_matrices                       = CR_matrices;
% 
%     results.(fieldNames{horizon_i}).cRMSEs                            = cRMSEs;
%     results.(fieldNames{horizon_i}).cMAEs                             = cMAEs;
%     results.(fieldNames{horizon_i}).cCRPS                             = cCRPS;
% 
%     results.(fieldNames{horizon_i}).GWtests_all_models_RMSE           = GWtests_all_models_RMSE;
%     results.(fieldNames{horizon_i}).GWtests_all_models_RMSE_09q1_15q4 = GWtests_all_models_RMSE_09q1_15q4;
%     results.(fieldNames{horizon_i}).GWtests_all_models_RMSE_16q1_19q4 = GWtests_all_models_RMSE_16q1_19q4;
% 
%     results.(fieldNames{horizon_i}).GWtests_all_models_MAE             = GWtests_all_models_MAE;
%     results.(fieldNames{horizon_i}).GWtests_all_models_MAE_09q1_15q4   = GWtests_all_models_MAE_09q1_15q4;
%     results.(fieldNames{horizon_i}).GWtests_all_models_MAE_16q1_19q4   = GWtests_all_models_MAE_16q1_19q4;
% 
%     results.(fieldNames{horizon_i}).GWtests_all_models_CRPS            = GWtests_all_models_CRPS;
%     results.(fieldNames{horizon_i}).GWtests_all_models_CRPS_09q1_15q4  = GWtests_all_models_CRPS_09q1_15q4;
%     results.(fieldNames{horizon_i}).GWtests_all_models_CRPS_16q1_19q4  = GWtests_all_models_CRPS_16q1_19q4;
% 
% end %for model_i
% 
% end
% 
