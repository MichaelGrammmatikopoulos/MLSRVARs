%======================================================================
%
%          MLSRVARs : MACHINE LEARNING SHADOW RATE VARs                        
%
%          Date : MARCH, 2025
%
%          Code Written By: Michael Grammatikopoulos    
%                    email: Michael.Grammatikopoulos@moodys.com  
%                           
%======================================================================
% 
% The working paper and supplementary appendices are available here:
% https://github.com/MichaelGrammatikopoulosMoodys
% 
% This code comes without technical support of any kind.  It is expected 
% to reproduce the results reported in the paper. Under no circumstances 
% will the author be held responsible for any use (or misuse) of this code 
% in any way.
%
% The views in this paper are solely of the author and do not represent 
% the views of Moody's Analytics or the Moody’s Corporation
%
%======================================================================
%
% Main citations: 
% 
% A.) Carriero, A., Clark T., Marcellino, M., Mertens, E.: 
% Shadow-Rate VARs (2024)
% Working Paper, Federal Reserve Bank of Cleveland
% https://www.elmarmertens.com/
% functions updated as of January 2024
%
% B.) Gefang, D., Koop, G. and Poon, A.:
% Forecasting Using Variational Bayesian Inference in Large Vector 
% Autoregressions with Hierarchical Shrinkage (2023), 
% International Journal of Forecasting, vol 39(1), pp. 346-363.
% https://sites.google.com/view/aubreybcpoon/research
%
% This code employs the Block Hybrid variant of the Shadow rate. The 
% equations for macroeconomic variables only include lagged actual rates 
% on the right-hand side. This is based on the premise that actual interest
% rates should influence spending and investment decisions. Conversely, the 
% presence of lagged shadow rates in the model's interest rate block aids 
% in encapsulating the lower-for-longer or make-up strategies in monetary 
% policy when predicting future interest rates.
%
% MLSRVARs Revision History:
% Februray 27, 2024 - Version 1
% September 1, 2024 - Version 2
% March 20, 2025 - Version 3
%======================================================================

%% ===================== HOUSEKEEPING =====================

close all
clearvars -global
clear; clc;

% Call the programs that loads raw data vintages
cd_path = 'YOURFOLDER/MLSVARs_main/';
pathdata = [cd_path, '/data'];
pathfunctions = [cd_path,'/matlab_codes/functions'];
pathfunctions2 = [cd_path,'/matlab_codes/functions/export_fig'];
addpath(pathdata);
addpath(pathfunctions);
addpath(pathfunctions2);
addpath('C:/Program Files/MATLAB/2023b/toolbox/econ/econ')
addpath('YOURFOLDER/MLSVARs_main/data\WuXia')

% Create SSP dataset
make_ssp_dataset

clearvars -except pathfunctions pathdata cd_path OHS_hyper_set idx spf_dataset_SSP

WUXIAshadow = matlab.io.datastore.FileSet([pathdata,'/WuXia/WuXiaShadowRate.xlsx']);
WUXIAshadow2 = spreadsheetDatastore(WUXIAshadow);
WUXIAshadow2.Sheets=1;
WUXIAshadow2 = read(WUXIAshadow2);

data_load

% Apply transformations and save the transformed vintages
data_transform

% initialising the random number generator
rng(2,"twister")

%% Global options
FLOOR_ELB                   = 0.25;                      % Raw non scaled value of FFR lower bound
CHOOSE_VAR                  = 4 ;

if  CHOOSE_VAR == 1
    which_VAR = 'small_VAR';
elseif  CHOOSE_VAR == 2
    which_VAR = 'financial_VAR';                         % 1: SMALL MACRO VAR: GDP, UNEMP, CPI, FEDFUNDS;
elseif  CHOOSE_VAR == 3                                  % 2: FINANCIAL VAR: the above plus:  6M ,  1Y , 5Y ,  10Y ,  BAA  ;
    which_VAR = 'macro_VAR';                             % 3: MEDIUM MACRO VAR: GDP, UNEMP, CPI, FEDFUNDS, INDPRO, MCUMFN , EXUSUK , M2SL , PINCOME , PCECC96 , PPIACO
elseif  CHOOSE_VAR == 4                                  % 4: FULL VAR
    which_VAR = 'full_VAR';
end

ROLLING_WINDOW              = 0;                         % 1: Rolling window; 0: Expanding window
NORMALIZE_DATA              = 0;                         % 1: Scale the Y ~ N(μ,σ^2) dataset to Z = (Y-μ)/σ ~ N(0,1)
nloop                       = 3000;                     % no of draws
burnin                      = nloop/3;                   % no of burn-in draws
elb_gibbs_burn              = 100;

p                           = 2;                         % no. of lags
forecast_horizon_set        = [1 4 8];
forecast_horizon            = max(forecast_horizon_set); % number of steps forecasted
forecast_Ndraws             = 1*(nloop-burnin);          % draws sampled from the predictive density
draw_predictive_densities   = 1;
thin_gibbs                  = 1;                         % draws for simulation.
check_stationarity          = 0;                         % truncatenonstationary draws?
shrink_all_to_zero          = 1;

IRF1scale                   = 1;
draw_GIRFs                  = 1;
irf_forecast_horizon        = 16;
variable_to_shock           = 'FEDFUNDS'; % CPIAUCSL

only_homoskedastic_models   = 0;
only_stvol_models           = 0;
only_srp_models             = 0;
only_ssp_models             = 0;
only_DL_models              = 0;
only_BLASSO_models          = 0;
only_SSVS_models            = 0;

prc70                       = normcdf([-1 1]) * 100;
Estrella_identification     = 0;
train_VARp_1_trainAR1_0     = 0;
VARp_L_prior                = 0;
only_pre_covid_samples      = 1;

%% All models loop preparation
% The program works without being case sensitive, and placement of the
% does not matter, i.e. 'stvol_ssvs' == 'SSVS_StVol'

% Some error messages for stability of the toolbox
if only_homoskedastic_models == 1 && only_stvol_models == 1
    error('Both only_homoskedastic_models and only_stvol_models are 1. The program will stop.')
end

if only_stvol_models+only_stvol_models+only_SSVS_models > 1 
    error('When you isolate models, only one of the shrinkage methods (DL, BLASSO or SSVS) can be selected')
end

%% Prior mnemonics:
% StVol:  Stochastic Volatility
% SSP:    Steady State Prior                   | SRP:    Shadow Rate Prior  
% SSVS:   Stochastic Search Variable Selection | A-BLASSO: Adaptive Bayesian Least Absolute Shrinkage & Selection Operator 

all_models = [
    "minnesota",     "minnesota_ssp",      "minnesota_stvol",      "minnesota_stvol_ssp",...
    "minnesota_srp", "minnesota_ssp_srp",  "minnesota_stvol_srp",  "minnesota_stvol_ssp_srp",...
    "ssvs",          "ssvs_ssp",           "ssvs_stvol",           "ssvs_stvol_ssp",...
    "ssvs_srp",      "ssvs_ssp_srp",       "ssvs_stvol_srp",       "ssvs_stvol_ssp_srp",...
    "blasso_A",      "blasso_A_ssp",       "blasso_A_stvol",       "blasso_A_stvol_ssp",...
    "blasso_A_srp",  "blasso_A_ssp_srp",   "blasso_A_stvol_srp",   "blasso_A_stvol_ssp_srp",    ...
    "DirLap",        "DirLap_ssp",         "DirLap_stvol",         "DirLap_stvol_ssp",...
    "DirLap_srp",    "DirLap_ssp_srp",     "DirLap_stvol_srp",     "DirLap_stvol_ssp_srp"
    ];

models_pretty

if only_stvol_models==1
    ndx_stvol_models=find(contains(all_models,'stvol'));
    all_models=all_models(ndx_stvol_models);
    all_models_pretty = all_models_pretty(ndx_stvol_models);
    all_models_backup = all_models;
end

if only_srp_models==1
    ndx_srp_models=find(contains(all_models,'srp'));
    all_models=all_models(ndx_srp_models);
    all_models_pretty = all_models_pretty(ndx_srp_models);
end

if only_ssp_models==1
    ndx_ssp_models=find(contains(all_models,'ssp'));
    all_models=all_models(ndx_ssp_models);
    all_models_pretty = all_models_pretty(ndx_ssp_models);
end

if only_DL_models==1
    ndx_DL_models=find(contains(all_models,'DirLap'));
    all_models=all_models(ndx_DL_models);
    all_models_pretty = all_models_pretty(ndx_DL_models);
end

if only_BLASSO_models==1
    ndx_BLASSO_models=find(contains(all_models,'blasso'));
    all_models=all_models(ndx_BLASSO_models);
    all_models_pretty = all_models_pretty(ndx_BLASSO_models);
end

if only_SSVS_models==1
    ndx_SSVS_models=find(contains(all_models,'ssvs'));
    all_models=all_models(ndx_SSVS_models);
    all_models_pretty = all_models_pretty(ndx_SSVS_models);
end

if only_homoskedastic_models==1
    ndx_homoskedastic_models=find(~contains(all_models,'stvol'));
    all_models=all_models(ndx_homoskedastic_models);
    all_models_pretty = all_models_pretty(ndx_homoskedastic_models);
    all_models_backup = all_models;
end

% In case you selected a subset of models save them here
all_models_subset = all_models;

% for simplicity, put first model as benchmark
benchmark_model = all_models(1);
ndxBENCH        = find(ismember(all_models, benchmark_model));
no_of_models    = size(all_models,2);

% Set the number of workers to 16
delete(gcp('nocreate'))
numWorkers = no_of_models;
parpool('local', numWorkers);

%% Set up the Forecast Evaluation (FE) samples
% For h = 8, forecasting 2009Q1 in 2007Q1, forecasting 2019Q4 in 2017Q4.
% For h = 4, forecasting 2009Q1 in 2008Q1, forecasting 2019Q4 in 2018Q4.
% For h = 1, forecasting 2009Q1 in 2008Q4, forecasting 2019Q4 in 2019Q3.

first_actual     = '2009 Q1';
last_actual      = '2019 Q4';
idx_first_actual = find(ismember(string(quarter_names_initial),{first_actual}));
idx_last_actual  = find(ismember(string(quarter_names_initial),{last_actual}));
first_vintage    = '2007 Q1';
last_vintage     = '2019 Q3';
idx_first_vintage = find(ismember(string(quarter_names_initial),{first_vintage}));
idx_last_vintage  = find(ismember(string(quarter_names_initial),{last_vintage}));
first_forecast   = '2007 Q2';
last_forecast    = '2019 Q4';
idx_first_forecast = idx_first_vintage+1;
idx_last_forecast  = idx_last_vintage+1;

quarter_names_FEsample        = quarter_names_initial(idx_first_actual:idx_last_actual,:);
quarters_forecast             = quarter_names_initial(idx_first_forecast:idx_last_forecast,:);
quarter_names_final           = quarter_names_initial(idx_first_vintage:idx_last_vintage);

all_forecast_samples = {};
all_forecast_samples_idx = {};
for ii = 1:length(forecast_horizon_set)
    forecast_horizon_i = forecast_horizon_set(ii);
    forecast_sample_i = [quarter_names_initial(idx_first_actual-forecast_horizon_i) quarter_names_initial(idx_last_actual-forecast_horizon_i)];
    all_forecast_samples{ii}=forecast_sample_i;
    all_forecast_samples_idx{ii}=find(ismember(string(quarters_forecast),string([quarter_names_initial(idx_first_actual-forecast_horizon_i+1) quarter_names_initial(idx_last_actual-forecast_horizon_i+1)])));
end

FEsample                            = find(ismember(string(quarter_names_initial),{first_vintage, last_vintage}));

% Adjust FEsample based on the forecast horizon (?) Not if you are running
% multiple h steps ahead at the same time. Example, if you run the code for
% h = 4, you have already produced h = 1 step ahead forecasts.

FEsample(2)                         = FEsample(2);
no_of_samples                       = FEsample(2)-FEsample(1)+1;
% ndxfcst_horizon                     = forecast_horizon:forecast_horizon:forecast_horizon*no_of_samples;

%% First Print Actuals

FPActuals_sample                    = find(ismember(string(quarter_names_initial),{first_actual, last_actual}));
FPActuals                           = cell2mat(Y_raw_Transformed_vintages(FPActuals_sample(1)));
FPActuals                           = FPActuals(end,:);

for i = 2:(FPActuals_sample(2)-FPActuals_sample(1)+1)
    FPActuals_i                     = cell2mat(Y_raw_Transformed_vintages(FPActuals_sample(1)+i-1));
    FPActuals_i                     = FPActuals_i(end,:);
    FPActuals                       = [FPActuals;FPActuals_i];
end

FEstruct = struct();
FEstruct.first_vintage=first_vintage;
FEstruct.last_vintage=last_vintage;
FEstruct.idx_first_vintage=idx_first_vintage;
FEstruct.idx_last_vintage=idx_last_vintage;
FEstruct.quarter_names_initial=quarter_names_initial;
FEstruct.FPActuals=FPActuals;

FPActuals_sample_irf                    = find(ismember(string(quarter_names_initial),{'2007 Q2', '2019 Q4'}));
FPActuals_irf                           = cell2mat(Y_raw_Transformed_vintages(FPActuals_sample_irf(1)));
FPActuals_irf                           = FPActuals_irf(end,:);

for i = 2:(FPActuals_sample_irf(2)-FPActuals_sample_irf(1)+1)
    FPActuals_i_irf                    = cell2mat(Y_raw_Transformed_vintages(FPActuals_sample_irf(1)+i-1));
    FPActuals_i_irf                    = FPActuals_i_irf(end,:);
    FPActuals_irf                       = [FPActuals_irf;FPActuals_i_irf];
end

%% Select Variable set
if CHOOSE_VAR == 1 % SMALL MACRO VAR:
    ndxMODEL = find(ismember(var_mnemonic,{'CPIAUCSL','GDP','UNRATE',...
        'FEDFUNDS'}));
elseif CHOOSE_VAR == 2 % FINANCIAL VAR
    ndxMODEL = find(ismember(var_mnemonic,{'CPIAUCSL','GDP','UNRATE',...
        'FEDFUNDS','TB3MS','TB6MS','GS1','GS3','GS5','GS10','BAA'}));
elseif CHOOSE_VAR == 3 % MEDIUM MACRO VAR
    ndxMODEL = find(ismember(var_mnemonic,{'CPIAUCSL','GDP','UNRATE',...
        'PPIACO','PINCOME','PAYEMS','INDPRO','HOUST','PCECC96','M2SL','MCUMFN','EXUSUK',...
        'FEDFUNDS'}));
elseif CHOOSE_VAR == 4 % FULL VAR
    ndxMODEL = find(ismember(var_mnemonic,{'CPIAUCSL','GDP','UNRATE',...
        'PPIACO','PINCOME','PAYEMS','INDPRO','HOUST','PCECC96','M2SL','MCUMFN','EXUSUK',...
        'FEDFUNDS','TB3MS','TB6MS','GS1','GS3','GS5','GS10','BAA','INFEXP'}));
end

%% Scale data?
dataT_original = Y_raw_Transformed_vintages{1,1};
dataT_original = dataT_original(3:end, ndxMODEL);

%% Create inflation expectations proxy by using the spread of the 10y and 1y bonds.
var_mnemonic_i = var_mnemonic(ndxMODEL);

if NORMALIZE_DATA == 1
    % Scale the original data without transformations
    mi_unscale = mean(dataT_original);
    sigma_unscale = std(dataT_original);
    ZdataT = (dataT_original-mean(dataT_original))./std(dataT_original);
else
    ZdataT = dataT_original;
end
Y0 = ZdataT(1:p,:);
Y = ZdataT(p+1:end,:);
plot(Y)
[T,N]=size(Y);
Klagreg    = N*p;
kappa_mikro = N+ p*N^2;
mi_mikro = N*(N-1)/2;

tcode = tcode_full(ndxMODEL);
pretty_names = fredMDprettylabel(var_mnemonic_i);

cumcode = logical(1*ones(N,1));
% Add inflation expectations to the cumulative graphs
% cumcode(~ismember(var_mnemonic_i,'INFEXP') & (tcode == 1 | tcode==4))=0;
cumcode((tcode == 1 | tcode==4))=0;

%% ==================== HYPERPARAMETERS ====================

hyperparameters_minnesota = [0.05 0.5 100 2];
L_const_min = 10;

%mean for Minnesota prior
setMinnesotaMean2 % 1 for variables in levels, 0 for all else
minnesotaPriorMean = NaN(N,1);
for n = 1 : N
    switch tcode(n)
        case {1,4}
            minnesotaPriorMean(n) = 1;
        otherwise
            minnesotaPriorMean(n) = 0;
    end
end

% Ensure SSP prior is 0 for the first and second differenced variables
SSP_id = ones(N,1);
for n = 1 : N
    switch tcode(n)
        case {2,3,6}
            SSP_id(n) = 0;
    end
end

% Prior on conditional mean coefficients, use Minnesota setup
if train_VARp_1_trainAR1_0 == 1
    % VAR(p) training
    VARresid=NaN(T-p,N);
    % for i=1:N
    tempY                       = ZdataT;
    Xt0                         = [];
    lags                        = p;
    for i=1:lags
        Xt0(:,(i-1)*N+1:i*N)    = tempY(lags-i+1:end-i,:);
    end
    Yt1 = ZdataT(lags+1:end,:);
    Xt0=[ones(size(Xt0,1),1) Xt0];

    VARresid=Yt1-Xt0*(Xt0\Yt1);
    % beta_prior_mean_training = (Xt0\Yt1);
    % minnesotaPriorMean = diag(beta_prior_mean_training(2:end,:));

    VAR_VCOV = VARresid'*VARresid / (T-size(Xt0,2));
    VAR_s2= diag(diag(VAR_VCOV));

    % [invA0prior Sigmaprior] = ldl(VAR_VCOV);
    invA0prior = inv(chol(VAR_VCOV,'lower')*inv(diag(diag(chol(VAR_VCOV,'lower')))));
    L_prior_means = tril(inv(invA0prior),-1)';

    thetas_L = cell(N-1,1);
    for j = 2:N
        thetas_L{j} =  0.1*L_prior_means(1:j-1,j);
        %     thetas_L{j} =  0*L_prior_means(1:j-1,j);
    end

    Pi_pm=zeros(N * Klagreg,1); Pi_pv=eye(N * Klagreg); co=0;
    % todo: Pi_pv could be encoded as vector, only using diagonal elements anyhow
    sigma_const = NaN(1,N);
    for i=1:N
        sigma_const(i)=VAR_s2(i,i)*hyperparameters_minnesota(3);  % this sets the prior variance on the intercept
        for l=1:p; %#ok<*NOSEL>
            for j=1:N
                co=co+1;
                if (i==j)
                    if l==1
                        Pi_pm(co)=minnesotaPriorMean(i); % this sets the prior means for the first own lag coefficients.
                    end
                    Pi_pv(co,co)=hyperparameters_minnesota(1)/(l^hyperparameters_minnesota(4)); % prior variance, own lags
                else
                    Pi_pv(co,co)=(VAR_s2(i,i)/VAR_s2(j,j)*hyperparameters_minnesota(1)*hyperparameters_minnesota(2)/(l^hyperparameters_minnesota(4))); % prior variance, cross-lags
                end
            end
        end
    end
else
    % AR(1) training
    beta_ij = [];
    for i=1:N
        yt_0=[ones(T-1,1) Y(1:end-1,i) ];
        yt_1=Y(2:end,i);
        ARresid(:,i)=yt_1-yt_0*(yt_0\yt_1);
    end

    AR_s2= diag(diag(ARresid'*ARresid))./(T-2);
    Pi_pm=zeros(N * Klagreg,1); Pi_pv=eye(N * Klagreg); co=0;
    % todo: Pi_pv could be encoded as vector, only using diagonal elements anyhow
    sigma_const = NaN(1,N);
    for i=1:N
        sigma_const(i)=AR_s2(i,i)*hyperparameters_minnesota(3);  % this sets the prior variance on the intercept
        for l=1:p; %#ok<*NOSEL>
            for j=1:N
                co=co+1;
                if (i==j)
                    if l==1
                        Pi_pm(co)=minnesotaPriorMean(i); % this sets the prior means for the first own lag coefficients.
                    end
                    Pi_pv(co,co)=hyperparameters_minnesota(1)/(l^hyperparameters_minnesota(4)); % prior variance, own lags
                else
                    Pi_pv(co,co)=(AR_s2(i,i)/AR_s2(j,j)*hyperparameters_minnesota(1)*hyperparameters_minnesota(2)/(l^hyperparameters_minnesota(4))); % prior variance, cross-lags
                end
            end
        end
    end
end

% Prior hyperparameters
nu = size(ndxMODEL,2); S0 = 0.01*nu; 
Vbeta = 1; Va= 1;

% SSVS without Minnesota scaling
tau0_c    = 0.0001;
tau1_c    = 1/9;

tau0_L    = tau0_c;
tau1_L    = tau1_c;

probability_i  = 1/2;

% SSVS with Minnesota
tau0_min = 0.5;
tau1_min = 2;

% Stochastic Volatility priors
nu0 = 5; S0h = 0.01;
ah  = 0; Vh  = 10;

% DL priors
abeta = 0.5; agam = 0.5; %abeta = 0.95; agam = 0.95;

% BLASSO
a0_c = N;                b0_c = 1;
a0_L = N;                b0_L = 1;

a0_global = N*(N-1)/2;  b0_global   = 1/N;

% SSP-SSVS
load VAR_s2.mat
tau0_MU = 0.0001*diag(VAR_s2(:,ndxMODEL))./diag(VAR_s2(:,ndxMODEL));
tau1_MU = 10;

ndx_variable_to_shock = find(ismember(var_mnemonic_i,{variable_to_shock}));

shadowrateTails_all_models_samples  = cell(no_of_models,1);
IRFPlus_draws_all_models_samples    = cell(no_of_models,1);
IRFMinus_draws_all_models_samples   = cell(no_of_models,1);
% yforecast_all_models_samples        = cell(no_of_models,1);
fcstYdraws_all_models_samples       = cell(no_of_models,1);
post_h_all_models_samples           = cell(no_of_models,1);
MU_all_models_samples               = cell(no_of_models,1);

failures_all_models_samples         = cell(no_of_models,1);

%% keep only the wanted vintages for the FE
Y_vintages_final                    = Y_raw_Transformed_vintages(FEsample(1):FEsample(2),1);
Y_vintages_final_ROLLING            = [];
naming                              = quarter_names_initial(FEsample(1):FEsample(2),1);

% Settings table
ssvs_tbl = table(tau1_c,tau0_c);
ssvs_min_tbl = table(tau1_min,tau0_min);

Min_tbl = table(hyperparameters_minnesota(1), hyperparameters_minnesota(2), hyperparameters_minnesota(3), hyperparameters_minnesota(4), L_const_min);
Min_tbl=renamevars(Min_tbl,["Var1","Var2","Var3","Var4","L_const_min"],["lambda1","lambda2","lambda3","lambda4","lambda5"]);
ablasso_tbl = table(a0_c,b0_c,a0_L,b0_L);
gblasso_tbl = table(a0_global,b0_global);
dilap_tbl = table(abeta,agam);
sv_tbl = table(nu0,S0h,ah,Vh);

display_settings_HSSRVAR
start_time_parfor = clock;

parfor model_i = 1:no_of_models

    try 

    time_start = datetime('now');

    % Set up the model / prior type string
    prior_type                      = all_models(model_i)

    % allocate memory / storage
    shadowrateTails_all_samples = cell(no_of_samples,1);
    yforecast_all_samples       = cell(no_of_samples,1);
    IRFPlus_tails_all_samples   = cell(no_of_samples,1);
    IRFMinus_tails_all_samples  = cell(no_of_samples,1);
    fcstYdraws_all_samples      = cell(no_of_samples,1);
    post_h_all_samples          = cell(no_of_samples,1);
    MU_all_samples              = cell(no_of_samples,1);

    Y_vintages_final_ROLLING    = cell(no_of_samples,1);
    failures_all_samples        = cell(no_of_samples,1);
    
    %% MCMC starts here
    randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
    disp(append('Starting MCMC.... ',all_models(model_i)));
    disp(' ' );
    start_time = clock;
    list_of_vintages = 1:no_of_samples;

    [MU_all_samples, post_h_all_samples, yforecast_all_samples, fcstYdraws_all_samples, IRFPlus_tails_all_samples, IRFMinus_tails_all_samples, shadowrateTails_all_samples] = parfor_vintages(prior_type, ...
    no_of_samples, list_of_vintages, Y_vintages_final, ndxMODEL, ...
    ROLLING_WINDOW, NORMALIZE_DATA, var_mnemonic_i, p, FLOOR_ELB, spf_dataset_SSP, forecast_Ndraws, ...
    nloop,burnin,forecast_horizon,irf_forecast_horizon,SSP_id,minnesotaPriorMean,sigma_const, ...
    Pi_pv,Pi_pm,VARp_L_prior,L_const_min,Estrella_identification,thin_gibbs,tau0_c, ...
    tau1_c,tau0_L,tau1_L,tau1_min,tau0_min,probability_i,nu0,S0h,ah,Vh,a0_c,b0_c,a0_L,b0_L,a0_global,b0_global, ...
    ndx_variable_to_shock,check_stationarity,draw_predictive_densities,prc70, ...
    forecast_horizon_set,shadowrateTails_all_samples, draw_GIRFs, IRF1scale, ...
    Klagreg, N, abeta, agam, elb_gibbs_burn, yforecast_all_samples, fcstYdraws_all_samples, ...
    IRFPlus_tails_all_samples, IRFMinus_tails_all_samples, post_h_all_samples, MU_all_samples, shrink_all_to_zero, ...
    tau0_MU, tau1_MU, S0, nu)

    shadowrateTails_all_models_samples{model_i,1}   = shadowrateTails_all_samples;
    IRFPlus_tails_all_models_samples{model_i,1}     = IRFPlus_tails_all_samples;
    IRFMinus_tails_all_models_samples{model_i,1}    = IRFMinus_tails_all_samples;
    fcstYdraws_all_models_samples{model_i,1}        = fcstYdraws_all_samples;
    for vint_i = 1:no_of_samples
        yforecast_all_models_samples{model_i,1}{vint_i,1}     = mean(fcstYdraws_all_samples{vint_i,1},3);
    end

    post_h_all_models_samples{model_i,1}            = post_h_all_samples;
    MU_all_models_samples{model_i,1}                = MU_all_samples;

    if contains(prior_type,'stvol','ignorecase',true)
        post_h_all_models_samples{model_i,1} = post_h_all_samples;
    end

    time_end = datetime('now');
    computational_time_all_models(model_i,1) = hours(time_end-time_start);
    disp( [append(all_models_pretty(model_i) ,': MCMC took ', num2str( etime( clock, start_time) ), ' seconds') ] )
    catch ME
        fprintf(ME.message);
    end

end %for model_i

actual_runtime = num2str( etime( clock, start_time_parfor) );
disp( ['All models` MCMCs took '  num2str( etime( clock, start_time_parfor) ) ' seconds' ] );

current_datetime = datetime('now');
formattedDate = string(datestr(current_datetime, 'dd_mmm_yyyy_HH_MM_SS'));
cd = 'T:\PROJECTS\_MG_docs\MLSVARs\revision1\';

metadata_dir                = append(cd,'output\',which_VAR,'\',formattedDate);
forecast_fancharts          = append(metadata_dir,'\forecast_fancharts\');
GIRFs_dir                   = append(metadata_dir,'\GIRFs\');
tables_dir                  = append(metadata_dir,'\tables\');
shadow_rates_graphs_dir     = append(metadata_dir,'\shadow_rates_graphs\');
ssp_graphs_dir              = append(metadata_dir,'\ssp_graphs\');
volatilities_dir            = append(metadata_dir,'\volatilities\');

mkdir(metadata_dir);
mkdir(forecast_fancharts);
mkdir(GIRFs_dir);
mkdir(shadow_rates_graphs_dir);
mkdir(ssp_graphs_dir);
mkdir(tables_dir);
mkdir(volatilities_dir);

save(append(metadata_dir,'\',which_VAR,'_',string(nloop),'_',regexprep(regexprep(string(datetime("now")), ' ', '_'),':','_'),'_results_prelim.mat'), '-v7.3');

% Create results
[results] = create_results(yforecast_all_models_samples, ...
    fcstYdraws_all_models_samples, FPActuals, N, no_of_models, ...
    ndxMODEL,no_of_samples, ndxBENCH, ...
    forecast_horizon_set,all_forecast_samples_idx,quarters_forecast);
make_strings_for_graphs

% Create tables and graphs

% IRF graph settings
vint1='2015 Q1';
vint2='2018 Q1';
vintages_to_run = [find(ismember(string(quarter_names_final),{vint1})) find(ismember(string(quarter_names_final),{vint2}))];
mnemonics_for_girfs = {'FEDFUNDS','HOUST','UNRATE','INFEXP', ...
                'INDPRO','PAYEMS','PPIACO','PCECC96','M2SL'};

cumul_IRFPlus_tails_all_models_samples = IRFPlus_tails_all_models_samples;
for sample_i = 1:no_of_samples
    for model_i = 1:no_of_models
        for quant_i = 1:3
            cumul_IRFPlus_tails_all_models_samples{model_i}{sample_i}(cumcode,:,quant_i) = cumsum(cumul_IRFPlus_tails_all_models_samples{model_i}{sample_i}(cumcode,:,quant_i),2);
        end
    end
end

cumul_IRFMinus_tails_all_models_samples = IRFMinus_tails_all_models_samples;
for sample_i = 1:no_of_samples
    for model_i = 1:no_of_models
        for quant_i = 1:3
            cumul_IRFMinus_tails_all_models_samples{model_i}{sample_i}(cumcode,:,quant_i) = cumsum(cumul_IRFMinus_tails_all_models_samples{model_i}{sample_i}(cumcode,:,quant_i),2);
        end
    end
end

[table_graph_results] = create_tables_and_graphs(naming,irf_forecast_horizon,no_of_models,draw_GIRFs, ndxBENCH,all_models_pretty, pretty_names, var_mnemonic_i, ...
    tcode, minnesotaPriorMean, which_VAR, results, ...
    tables_dir, GIRFs_dir, shadow_rates_graphs_dir, volatilities_dir, forecast_fancharts, CHOOSE_VAR,...
    IRFPlus_tails_all_models_samples, IRFMinus_tails_all_models_samples, vintages_to_run, mnemonics_for_girfs, ...
    all_models, Yraw_table_last_vintage, shadowrateTails_all_models_samples, ...
    FPActuals_irf, ndxMODEL, fcstYdraws_all_models_samples, ...
    post_h_all_models_samples, no_of_samples, spf_dataset_SSP, SSP_id, MU_all_models_samples,forecast_horizon_set,WUXIAshadow2);
