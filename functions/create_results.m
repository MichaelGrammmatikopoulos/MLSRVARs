function  [results] = create_results(yforecast_all_models_samples, ...
    fcstYdraws_all_models_samples, FPActuals, N, no_of_models, ndxMODEL, ...
    no_of_samples, ndxBENCH, forecast_horizon_set, all_forecast_samples_idx,quarters_forecast,window_i)

idx2015q4 = find(ismember(string(quarters_forecast),{'2015 Q4'}));
idx2019q4 = find(ismember(string(quarters_forecast),{'2019 Q4'}));

idxFE2015q4 = find(ismember(string(quarters_forecast),{'2015 Q4'}));

results = struct();
fieldNames = {'h'+string(forecast_horizon_set(1)), 'h'+string(forecast_horizon_set(2)), 'h'+string(forecast_horizon_set(3))};
for horizon_i = 1:length(forecast_horizon_set)

    forecast_horizon_i = forecast_horizon_set(horizon_i);
    range_i = all_forecast_samples_idx{horizon_i}(1):all_forecast_samples_idx{horizon_i}(2);

    rRMSEs              = NaN(1,N,no_of_models);
    rRMSEs_09q1_15q4    = NaN(1,N,no_of_models);
    rRMSEs_16q1_19q4    = NaN(1,N,no_of_models);

    rMAEs               = NaN(1,N,no_of_models);
    rMAEs_09q1_15q4     = NaN(1,N,no_of_models);
    rMAEs_16q1_19q4     = NaN(1,N,no_of_models);

    % Collect
    RMSE_matrices                           = NaN(1,N,no_of_models);
    RMSE_matrices_09q1_15q4                 = NaN(1,N,no_of_models);
    RMSE_matrices_16q1_19q4                 = NaN(1,N,no_of_models);
    MAE_matrices                            = NaN(1,N,no_of_models);
    MAE_matrices_09q1_15q4                  = NaN(1,N,no_of_models);
    MAE_matrices_16q1_19q4                  = NaN(1,N,no_of_models);
    full_forecast_matrices                  = {};

    % Collect cumulatives: cRMSEs, cMAEs
    cRMSEs              = NaN(size(FPActuals(:,ndxMODEL),1),size(FPActuals(:,ndxMODEL),2),no_of_models);
    cMAEs               = NaN(size(FPActuals(:,ndxMODEL),1),size(FPActuals(:,ndxMODEL),2),no_of_models);
    cCRPS               = NaN(size(FPActuals(:,ndxMODEL),1),size(FPActuals(:,ndxMODEL),2),no_of_models);

    % Collect ratios
    rCRPS               = NaN(1,N,no_of_models);
    rCRPS_09q1_15q4     = NaN(1,N,no_of_models);
    rCRPS_16q1_19q4     = NaN(1,N,no_of_models);

    for model_i = 1:no_of_models

        full_forecast_matrix = cell2mat(yforecast_all_models_samples{model_i,1});
        full_forecast_matrix = full_forecast_matrix(:,horizon_i);
        full_forecast_matrix = reshape(full_forecast_matrix(:,1,:),N,[])';

        full_forecast_matrix = full_forecast_matrix(range_i,:);
        full_forecast_matrices{model_i} = full_forecast_matrix;

        % Quantities for the GW test
        SE_matrix  = (FPActuals(:,ndxMODEL)-full_forecast_matrix).^2;
        AE_matrix   = abs(FPActuals(:,ndxMODEL)-full_forecast_matrix);

        % Quantities for the ratios
        RMSE_matrix            = sqrt(mean((FPActuals(:,ndxMODEL)-full_forecast_matrix).^2));
        RMSE_matrix_09q1_15q4  = sqrt(mean((FPActuals(1:idx2015q4,ndxMODEL)-full_forecast_matrix(1:idx2015q4,:)).^2));
        RMSE_matrix_16q1_19q4  = sqrt(mean((FPActuals(idx2015q4+1:end,ndxMODEL)-full_forecast_matrix(idx2015q4+1:end,:)).^2));

        MAE_matrix             = mean(abs(FPActuals(:,ndxMODEL)-full_forecast_matrix));
        MAE_matrix_09q1_15q4   = mean(abs(FPActuals(1:idx2015q4,ndxMODEL)-full_forecast_matrix(1:idx2015q4,:)));
        MAE_matrix_16q1_19q4   = mean(abs(FPActuals(idx2015q4+1:end,ndxMODEL)-full_forecast_matrix(idx2015q4+1:end,:)));

        % save all quantities
        RMSE_matrices(:,:,model_i)           = RMSE_matrix;
        RMSE_matrices_09q1_15q4(:,:,model_i) = RMSE_matrix_09q1_15q4;
        RMSE_matrices_16q1_19q4(:,:,model_i) = RMSE_matrix_16q1_19q4;

        MAE_matrices(:,:,model_i)            = MAE_matrix;
        MAE_matrices_09q1_15q4(:,:,model_i)  = MAE_matrix_09q1_15q4;
        MAE_matrices_16q1_19q4(:,:,model_i)  = MAE_matrix_16q1_19q4;

        SE_matrices(:,:,model_i) = SE_matrix;
        AE_matrices(:,:,model_i) = AE_matrix;

        fcstYdraws_all_samples          = cell(no_of_samples,1);
        fcstYdraws_all_models_samples_i = fcstYdraws_all_models_samples{model_i,1};
        for sample_i = 1:no_of_samples
            fcstYdraws_all_samples{sample_i,1} = fcstYdraws_all_models_samples_i{sample_i,1};
        end

        rRMSEs(:,:,model_i)             = round(RMSE_matrices(:,:,model_i)./RMSE_matrices(:,:,ndxBENCH),2);
        rRMSEs_09q1_15q4(:,:,model_i)   = round(RMSE_matrices_09q1_15q4(:,:,model_i)./RMSE_matrices_09q1_15q4(:,:,ndxBENCH),2);
        rRMSEs_16q1_19q4(:,:,model_i)   = round(RMSE_matrices_16q1_19q4(:,:,model_i)./RMSE_matrices_16q1_19q4(:,:,ndxBENCH),2);

        rMAEs(:,:,model_i)              = round(MAE_matrices(:,:,model_i)./MAE_matrices(:,:,ndxBENCH),2);
        rMAEs_09q1_15q4(:,:,model_i)    = round(MAE_matrices_09q1_15q4(:,:,model_i)./MAE_matrices_09q1_15q4(:,:,ndxBENCH),2);
        rMAEs_16q1_19q4(:,:,model_i)    = round(MAE_matrices_16q1_19q4(:,:,model_i)./MAE_matrices_16q1_19q4(:,:,ndxBENCH),2);

        %% CRPS
        actuals_i=0;
        for vint_i=range_i
            actuals_i=actuals_i+1;
            for nnn = 1 : N % loop over elements of Y
                CRPSscore(nnn,actuals_i) = crpsDraws(FPActuals(actuals_i,nnn), reshape(fcstYdraws_all_models_samples{model_i}{vint_i}(nnn,horizon_i,:),[],1),[]);
            end
        end
        CRPS_matrices(:,:,model_i) = CRPSscore';
        CRPS_matrices_09q1_15q4(:,:,model_i) = CRPSscore(:,1:idxFE2015q4)';
        CRPS_matrices_16q1_19q4(:,:,model_i) = CRPSscore(:,idxFE2015q4+1:end)';

        cRMSEs(:,:,model_i)   = cumsum(SE_matrices(:,:,model_i));
        cMAEs(:,:,model_i)    = cumsum(AE_matrices(:,:,model_i));
        cCRPS(:,:,model_i)    = cumsum(CRPS_matrices(:,:,model_i));

        rCRPS(:,:,model_i)              = round(CRPS_matrices(end,:,model_i)./CRPS_matrices(end,:,ndxBENCH),2);
        rCRPS_09q1_15q4(:,:,model_i)    = round(CRPS_matrices_09q1_15q4(end,:,model_i)./CRPS_matrices_09q1_15q4(end,:,ndxBENCH),2);
        rCRPS_16q1_19q4(:,:,model_i)    = round(CRPS_matrices_16q1_19q4(end,:,model_i)./CRPS_matrices_16q1_19q4(end,:,ndxBENCH),2);


        %% Perform the Giacomini-White test
        if model_i~=ndxBENCH
            for col_i = 1:N
                GWtests_all_models_RMSE(model_i,col_i)           = CPAtest(SE_matrices(:,col_i,model_i),SE_matrices(:,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
                GWtests_all_models_RMSE_09q1_15q4(model_i,col_i) = CPAtest(SE_matrices(1:idx2015q4,col_i,model_i),SE_matrices(1:idx2015q4,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
                GWtests_all_models_RMSE_16q1_19q4(model_i,col_i) = CPAtest(SE_matrices(idx2015q4+1:end,col_i,model_i),SE_matrices(idx2015q4+1:end,col_i,ndxBENCH),forecast_horizon_i,0.05,1);

                GWtests_all_models_MAE(model_i,col_i)            = CPAtest(AE_matrices(:,col_i,model_i),AE_matrices(:,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
                GWtests_all_models_MAE_09q1_15q4(model_i,col_i)  = CPAtest(AE_matrices(1:idx2015q4,col_i,model_i),AE_matrices(1:idx2015q4,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
                GWtests_all_models_MAE_16q1_19q4(model_i,col_i)  = CPAtest(AE_matrices(idx2015q4+1:end,col_i,model_i),AE_matrices(idx2015q4+1:end,col_i,ndxBENCH),forecast_horizon_i,0.05,1);

                GWtests_all_models_CRPS(model_i,col_i)           = CPAtest(cCRPS(:,col_i,model_i), cCRPS(:,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
                GWtests_all_models_CRPS_09q1_15q4(model_i,col_i) = CPAtest(cCRPS(1:idx2015q4,col_i,model_i), cCRPS(1:idx2015q4,col_i,ndxBENCH),forecast_horizon_i,0.05,1);
                GWtests_all_models_CRPS_16q1_19q4(model_i,col_i) = CPAtest(cCRPS(idx2015q4+1:end,col_i,model_i), cCRPS(idx2015q4+1:end,col_i,ndxBENCH),forecast_horizon_i,0.05,1);

            end
        end

    end %for model_i

    results.(fieldNames{horizon_i}).rRMSEs_09q1_15q4                  = reshape(rRMSEs_09q1_15q4,N,no_of_models)';
    results.(fieldNames{horizon_i}).rRMSEs_16q1_19q4                  = reshape(rRMSEs_16q1_19q4,N,no_of_models)';
    results.(fieldNames{horizon_i}).rRMSEs                            = reshape(rRMSEs,N,no_of_models)';
    results.(fieldNames{horizon_i}).RMSE_matrices                     = RMSE_matrices;

    results.(fieldNames{horizon_i}).rMAEs_09q1_15q4                   = reshape(rMAEs_09q1_15q4,N,no_of_models)';
    results.(fieldNames{horizon_i}).rMAEs_16q1_19q4                   = reshape(rMAEs_16q1_19q4,N,no_of_models)';
    results.(fieldNames{horizon_i}).rMAEs                             = reshape(rMAEs,N,no_of_models)';
    results.(fieldNames{horizon_i}).MAE_matrices                      = MAE_matrices;

    results.(fieldNames{horizon_i}).cRMSEs                            = cRMSEs;
    results.(fieldNames{horizon_i}).cMAEs                             = cMAEs;
    results.(fieldNames{horizon_i}).cCRPS                             = cCRPS;

    results.(fieldNames{horizon_i}).rCRPS                             = reshape(rCRPS,N,no_of_models)';
    results.(fieldNames{horizon_i}).rCRPS_09q1_15q4                   = reshape(rCRPS_09q1_15q4,N,no_of_models)';
    results.(fieldNames{horizon_i}).rCRPS_16q1_19q4                   = reshape(rCRPS_16q1_19q4,N,no_of_models)';

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

% Recusrive means calculation
for horizon_i = 1:length(forecast_horizon_set)
    forecast_horizon_i = forecast_horizon_set(horizon_i);
    for model_i = 1:no_of_models
        results.(fieldNames{horizon_i}).Recursive_Means_RMSE(:,:,model_i) = (movmean(results.(fieldNames{horizon_i}).cRMSEs(:,:,model_i),window_i,'endpoints','discard')./movmean(results.(fieldNames{horizon_i}).cRMSEs(:,:,ndxBENCH),window_i,'endpoints','discard')-1)*100;
        results.(fieldNames{horizon_i}).Recursive_Means_MAE(:,:,model_i) = (movmean(results.(fieldNames{horizon_i}).cMAEs(:,:,model_i),window_i,'endpoints','discard')./movmean(results.(fieldNames{horizon_i}).cMAEs(:,:,ndxBENCH),window_i,'endpoints','discard')-1)*100;
        results.(fieldNames{horizon_i}).Recursive_Means_CRPS(:,:,model_i) = (movmean(results.(fieldNames{horizon_i}).cCRPS(:,:,model_i),window_i,'endpoints','discard')./movmean(results.(fieldNames{horizon_i}).cCRPS(:,:,ndxBENCH),window_i,'endpoints','discard')-1)*100;
    end
end
