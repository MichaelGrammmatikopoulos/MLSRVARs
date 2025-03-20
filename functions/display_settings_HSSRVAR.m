
disp('                                              ');
disp('                                              ');
disp('    MACHINE LEARNING SHADOW RATE VARs         ');
disp('                                              ');
disp('                                              ');

all_settings = struct();
all_settings.FLOOR_ELB=FLOOR_ELB;
all_settings.CHOOSE_VAR=CHOOSE_VAR;
all_settings.ROLLING_WINDOW=ROLLING_WINDOW;
all_settings.NORMALIZE_DATA=NORMALIZE_DATA;
all_settings.nloop=nloop;
all_settings.burnin=burnin;
all_settings.elb_gibbs_burn=elb_gibbs_burn;
all_settings.p=p;
all_settings.forecast_horizon_set=forecast_horizon_set;
all_settings.forecast_horizon=forecast_horizon;
all_settings.forecast_Ndraws=forecast_Ndraws;
all_settings.draw_predictive_densities=draw_predictive_densities;
all_settings.thin_gibbs=thin_gibbs;
all_settings.check_stationarity=check_stationarity;
all_settings.shrink_all_to_zero=shrink_all_to_zero;
all_settings.IRF1scale=IRF1scale;
all_settings.draw_GIRFs=draw_GIRFs;
all_settings.irf_forecast_horizon=irf_forecast_horizon;
all_settings.variable_to_shock=variable_to_shock;
all_settings.only_stvol_models=only_stvol_models;
all_settings.only_ssp_models=only_ssp_models;
all_settings.only_DL_models=only_DL_models;
all_settings.only_BLASSO_models=only_BLASSO_models;
all_settings.prc70=prc70;
all_settings.VARp_L_prior=VARp_L_prior;
all_settings.train_VARp_1_trainAR1_0=train_VARp_1_trainAR1_0;
all_settings.Estrella_identification=Estrella_identification;        
all_settings.benchmark_model=benchmark_model;
all_settings.var_mnemonic_i = var_mnemonic_i;
all_settings.Min_tbl=Min_tbl;
all_settings.ssvs_tbl=ssvs_tbl;
all_settings.ablasso_tbl=ablasso_tbl;
all_settings.dilap_tbl=dilap_tbl;
all_settings.sv_tbl=sv_tbl;
all_settings.tau0_MU = tau0_MU;
all_settings.tau1_MU = tau1_MU;

disp('    SETTINGS:         ');
disp('                                              ');
% disp(struct2table(all_settings, 'AsArray', true))
disp(all_settings)