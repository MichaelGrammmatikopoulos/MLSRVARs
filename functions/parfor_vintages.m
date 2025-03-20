function  [MU_all_samples, post_h_all_samples, yforecast_all_samples, fcstYdraws_all_samples, IRFPlus_tails_all_samples, IRFMinus_tails_all_samples, shadowrateTails_all_samples] = parfor_vintages(prior_type, ...
    no_of_samples, list_of_vintages, Y_vintages_final, ndxMODEL, ...
    ROLLING_WINDOW, NORMALIZE_DATA, var_mnemonic_i, p, FLOOR_ELB, spf_dataset_SSP, forecast_Ndraws, ...
    nloop,burnin,forecast_horizon,irf_forecast_horizon,SSP_id,minnesotaPriorMean,sigma_const, ...
    Pi_pv,Pi_pm,VARp_L_prior,L_const_min,Estrella_identification,thin_gibbs,tau0_c, ...
    tau1_c,tau0_L,tau1_L,tau1_min,tau0_min,probability_i,nu0,S0h,ah,Vh,a0_c,b0_c,a0_L,b0_L,a0_global,b0_global, ...
    ndx_variable_to_shock,check_stationarity,draw_predictive_densities,prc70, ...
    forecast_horizon_set,shadowrateTails_all_samples, draw_GIRFs, IRF1scale, ...
    Klagreg, N,  abeta, agam, elb_gibbs_burn, yforecast_all_samples, fcstYdraws_all_samples, ...
    IRFPlus_tails_all_samples, IRFMinus_tails_all_samples, post_h_all_samples, MU_all_samples, ...
    shrink_all_to_zero,tau0_MU, tau1_MU, S0, nu)

for vint_ii = 1:no_of_samples

    vint_i = list_of_vintages(vint_ii);

    if vint_i == 24
        vint_i;
    end

    data = Y_vintages_final{vint_i,1};
%     disp(['last historical: ' quarter_names_final(vint_i)])
    data=data((p+1):end,ndxMODEL); % remove na rows stemming from the first and second differences

    % remove initial observationss in each vintage for rolling window
    % (fixed # of observations)
    if ROLLING_WINDOW == 1
        data = data(vint_i:end,:);
        Y_vintages_final_ROLLING{vint_i,1}=data;
    end

    % subset the full arrays to the model arrays
    %     var_mnemonic_i = var_mnemonic(ndxMODEL);
    N = size(data,2); % no. of variables
    Y0 = data(1:p,:);
    Y = data(p+1:end,:);
    Kbvar                    = N * p + 1;  % number of regressors per equation
    K                        = Kbvar;
    eyeN                     = eye(N);     % save an identity matrix for quicker calculation of inverse
    T                        = size(Y,1);  % no. of observations
    k                        = N+ p*N^2;
    np                       = N*p;
    m                        = N*(N - 1)/2;
    X                        = zeros(T,N*p);
    tempY                    = [Y0; Y];
    for i=1:p
        X(:,(i-1)*N+1:i*N)   = tempY(p-i+1:end-i,:);
    end
    X                        = [ones(T,1) X];

    % update pointers
    [T,K]=size(X);
    Klagreg    = K - 1; % number of lag regressors (without intercept)
    ndxKlagreg = 1 + (1 : Klagreg); % location of the Klag regressors in X
    % store original data matrices
    % these are scaled or unscaled depending on NORMALIZED_DATA
    Y_original = Y;
    X_original = X;
    dataT_original = data;

    % DL priors
    DL_scaler = 0.1; % default is 0.1
    psibeta = DL_scaler*ones(k,1);
    zetabeta = DL_scaler;
    varthetabeta = DL_scaler*ones(k,1);

    psia = DL_scaler*ones(m,1);
    zetaa = DL_scaler;
    varthetaa = DL_scaler*ones(m,1);

    idn = [];
    for i = 1:N
        if i > 1
            idn = [idn;ones(i-1,1)*i]  ;
        end
    end

    % Shadow rate prior perliminaries
    ndxSHADOWRATE   = find(ismember(var_mnemonic_i, {'FEDFUNDS', 'TB3MS', 'TB6MS', 'GS1','GS3'}));
    % define index of vars that need to obey ELB (at least out of sample)
    ndxOTHERYIELDS  = find(ismember(var_mnemonic_i, {'GS5', 'GS10','BAA'}));
    ndxYIELDS       = union(ndxSHADOWRATE, ndxOTHERYIELDS);
    ndxYIELDLAGS    = cat(2, false, repmat(ismember(1:N, ndxYIELDS), 1, p));
    Nyields         = length(ndxYIELDS);

    % Scale data
    if NORMALIZE_DATA == 1
        % Scale the original data without transformations
        mi_unscale = mean(dataT_original);
        sigma_unscale = std(dataT_original);
        ZdataT = (dataT_original-mean(dataT_original))./std(dataT_original);

        data = ZdataT;

        % Scale ELBbound
        ELBbound = (FLOOR_ELB-mi_unscale(ndxSHADOWRATE))./sigma_unscale(ndxSHADOWRATE);
        ELBboundyields = ((FLOOR_ELB)-mi_unscale(ndxYIELDS))./sigma_unscale(ndxYIELDS);
    else
        ZdataT = dataT_original;
        data = ZdataT;

        ELBbound = 0.25;
        ELBboundyields = repmat(FLOOR_ELB,size(ndxYIELDS));
    end

    Y0 = data(1:p,:);
    Y = data(p+1:end,:);
    X                        = zeros(T,N*p);
    tempY                    = [Y0; Y];
    for i=1:p
        X(:,(i-1)*N+1:i*N)   = tempY(p-i+1:end-i,:);
    end

    % Add constant only for Non-SSP models
    if ~contains(lower(prior_type),'ssp')
        X                     = [ones(T,1) X];
    end

    % SSP prior initialization
    if contains(lower(prior_type),'ssp')
        bols=[ones(T,1) X]\Y; % needs to be calculated with the constant
        no_of_regression_coef = size(bols,1)-1;
    else
        bols=X\Y; % needs to be calculated with the constant
        no_of_regression_coef = size(bols,1);
    end

    eyeN = eye(N);
    if contains(lower(prior_type),'ssp')
        M0 = table2array(spf_dataset_SSP(vint_i,1+ndxMODEL))'.*SSP_id;
        V0 = 1e-6*eye(N);
        invV0 = inv(V0);

        %% SSPssvs
        SSPssvs = 1;
        probability_i_MU=0.5;

        if NORMALIZE_DATA == 1
            M0 = (M0-mi_unscale')./sigma_unscale';
        end
    end

    F = [bols(2:N*p+1,:)' ; eye(N*(p-1),N*p)]; %companion form
    eyeF = eye(size(F,2));
    C=zeros(size(F,1),1);
    C(1:N)=bols(1,:)';
    MU = repmat(mean(Y),1,2)';

    % Create data matrices
    [Nobs,N]                 = size(data);
    Nshadowrates             = length(ndxSHADOWRATE);
    Ydata                    = data;
    shadowrateData           = data(:,ndxSHADOWRATE);
    ndxELB                   = shadowrateData <= ELBbound;
    Xactual                  = X;
    if contains(lower(prior_type),'ssp')
        Xactual = [ones(T,1) Xactual];
    end
    ELBdummy = data(:,ndxSHADOWRATE) <= ELBbound;  % find all data that are below or equal the ELB
    startELB = find(any(ELBdummy,2), 1);
    if size(startELB,1) == 0
        elbT0=inf;
    else
        elbT0    = startELB - 1 - p;
    end

    actualrateBlock = ~ismember(1:N, ndxYIELDS);

    % generate state vector for forecast jumpoff
    Nstates       = K + Nyields * p; % track lagged actual rates as well
    Xjumpoff      = zeros(Nstates,1);
    Xjumpoff(1)   = 1;

    for l=1:p
        Xjumpoff(1+(l-1)*N+(1:N)) = data(Nobs-(l-1),1:N);
    end

    if contains(prior_type,'srp','ignorecase',true)
        % add lagged actual rates to jumpoff
        for l=1:p
            Xjumpoff(K+(l-1)*Nyields+(1:Nyields)) = data(Nobs-(l-1),Nyields);
        end

        ndxSHADOWRATELAGS = cat(2, false, repmat(ismember(1:N, ndxSHADOWRATE), 1, p)); % preprend by false for CONST
        ndxSHADOWRATELAGS_initial=ndxSHADOWRATELAGS;
    end

    %% allocate memory for out-of-sample forecasts
    Ndraws      = forecast_Ndraws / (nloop-burnin);

    if Ndraws~=1
        if mod(forecast_Ndraws, nloop) ~= 0
            error('Ndraws must be multiple of nloop')
        end
    end
    fcstYdraws = NaN(N,forecast_horizon, Ndraws, nloop-burnin); % see reshape near end of script
    fcstYdraws_irf = NaN(N,irf_forecast_horizon, Ndraws, nloop-burnin); % see reshape near end of script
    fcstYdraws_irf_graphs = NaN(N,irf_forecast_horizon, Ndraws, nloop-burnin); % see reshape near end of script
    [fcstYdraws1plus, fcstYdraws1minus] = deal(NaN(N,forecast_horizon, Ndraws, nloop-burnin));

    %% prepare state space for forecasting

    if contains(prior_type,'srp','ignorecase',true)
        fcstA                  = zeros(Nstates,Nstates);
        fcstA(1,1)             = 1; % unit root for the constant
        fcstA(1+N+1:K,2:K)     = [eye(N*(p-1)),zeros(N*(p-1),N)]; % fill in lower part of companion form
        fcstA(K+Nyields+(1:Nyields*(p-1)),K+(1:Nyields*(p-1))) = eye(Nyields*(p-1)); % companion for actual rates
        ndxfcstActual          = K+(1:Nyields);
        ndxfcstShadow          = 1+ndxYIELDS;

        ndxfcstY               = 1+(1:N);
        fcstB                  = zeros(Nstates,N);
        fcstB(ndxfcstY,:)      = eye(N);

        %% EM, new: prepare ELB data
        shadowrateData(ndxELB)   = NaN;
        Ydata(:,ndxSHADOWRATE)   = shadowrateData;

        % yNaN
        YdataNaN                 = isnan(Ydata);
        Ydata(YdataNaN)          = 0;
        yNaN                     = YdataNaN(p+1:end,:);

        hasELBdata = T > elbT0;
        elbT       = max(0,T - elbT0);

        if hasELBdata
            elb.Nproposals     = 1e3;

            elb.gibbsburn      = elb_gibbs_burn;
            elb.ndxS           = ismember(1:N, ndxSHADOWRATE);

            elb.ndxY           = 1+(1:N);
            elb.T0             = elbT0;
            %             elb.T              = T - elb.T0;
            elb.yNaN           = yNaN(elb.T0+1:end,:)'; % note the transpose, needed for call of sampler later
            elb.sNaN           = elb.yNaN(elb.ndxS,:);
            if any(any(yNaN(1:elb.T0,ndxSHADOWRATE)))
                error('something off about elbT0')
            end

            if ~any(any(yNaN(elb.T0+1,ndxSHADOWRATE)))
                warning('elb.T0 + 1 should be a missing obs, but it is not ...')
            end
            elb.bound = ELBbound;
            % init transition matrix
            elb.A                = zeros(K,K);
            elb.A(1,1)           = 1; % unit root for the constant
            elb.A(1+N+1:end,2:end) = [eye(N*(p-1)),zeros(N*(p-1),N)]; % fill in lower part of companion form
            elb.A_backup = elb.A;
            elb.B            = zeros(K,N);
            % encode "prior" over initial conditions (which are fixed)
            if contains(prior_type,'ssp','ignorecase',true)
                elb.X0  = [1 X(elb.T0+1,:)]'; % Ensure a constant is included if you have SSP prior
            else
                elb.X0  = X(elb.T0+1,:)'; % time zero values of state. Recall that X already contains lagged values
            end
            % construct state vector for ELB state space
            dummy           = Ydata;
            dummy(YdataNaN) = NaN;
            dummy(:,ndxSHADOWRATE) = NaN; % to ignore all information about a actual or shadow rate
            elb.X   = ones(1+N*p,elbT);
            for l=0:p-1
                elb.X(1+N*l+(1:N),:) = dummy(p+elb.T0-l+(1:elbT),1:N)'; % note: T0 indexes into data after cutting out p lags
            end
        end
        tau = cell(N,1);
    else
        fcstA                  = zeros(K,K);
        fcstA(1,1)             = 1; % unit root for the constant
        fcstA(1+N+1:K,2:K)     = [eye(N*(p-1)),zeros(N*(p-1),N)]; % fill in lower part of companion form
        ndxfcstY               = 1+(1:N);
        fcstB                  = zeros(K,N);
        fcstB(ndxfcstY,:)      = eye(N);
    end

    % after creating all necessary matrices, recalculate quantities
    % without the constant term inside for the Steady State Prior
    if contains(prior_type,'ssp','ignorecase',true)
        K = N * p;
        k = p*N^2;

        % DL
        zetabeta = DL_scaler;
        varthetabeta = DL_scaler*ones(k,1);

    end

    if contains(prior_type,'ssp','ignorecase',true)
        OMEGA_pai   = sparse(diag(vec([reshape(diag(Pi_pv),Klagreg,N)]))); % prior variance of beta_store_mat
        MU_pai      = [reshape(Pi_pm,Klagreg,N)];                   % prior mean of beta_store_mat
    else
        OMEGA_pai   = sparse(diag(vec([sigma_const;reshape(diag(Pi_pv),Klagreg,N)]))); % prior variance of beta_store_mat
        MU_pai      = [zeros(1,N);reshape(Pi_pm,Klagreg,N)];                   % prior mean of beta_store_mat
    end
    iV   = sparse(diag(1./diag(OMEGA_pai)));
    invVthetas = reshape(getmaindiagonal(iV),[],N);

    theta0=cell(N,1);
    for iv = 1:N
        if iv ==1
            km = k/N;
        else
            km = k/N + iv - 1 ;
        end

        if VARp_L_prior == 1
            theta0{iv,1} = [zeros(k/N,1);thetas_L{iv}];
        else
            theta0{iv,1} = zeros(km,1);
        end

        theta0{iv,1}(1:K) = columnextract(MU_pai,iv);

        invVtheta{iv,1} = L_const_min*ones(km,1);
        invVtheta{iv,1}(1:K) = columnextract(invVthetas,iv);
        invVtheta{iv,1} = sparse(diag(invVtheta{iv,1}));

        invVtheta_min{iv,1} = sparse(diag(invVtheta{iv,1}));

        % Pure SSVS
        if contains(prior_type,'ssvs','ignorecase',true) && ~contains(prior_type,'_Min','ignorecase',true)
            if shrink_all_to_zero == 1
                theta0{iv,1} = 1e-6*ones(km,1);
            end

            invVtheta{iv,1} = 100*ones(km,1);;
            invVtheta{iv,1}(1:K) = 10;
            invVtheta{iv,1} = sparse(diag(invVtheta{iv,1}));
        end

        % DL
        if contains(prior_type,'DirLap','ignorecase',true)
            if shrink_all_to_zero == 1
                theta0{iv,1} = 1e-6*ones(km,1);
            end

            invVtheta{iv,1} = 100*ones(km,1);
            invVtheta{iv,1}(1:K) = 10;
            invVtheta{iv,1} = sparse(diag(invVtheta{iv,1}));
        end

        %BLASSO
        if contains(prior_type,'blasso','ignorecase',true)
            if shrink_all_to_zero == 1
                theta0{iv,1} = 1e-6*ones(km,1);
            end

            tau{iv,1} = 10^-6*ones(km,1);
            tau_c{iv} = tau{iv}(1:no_of_regression_coef);
            tau_L{iv} = tau{iv}(no_of_regression_coef+1:end);
            invVtheta{iv,1} = sparse(diag(1./tau{iv,1}));
            lambda{iv,1} = 1;
        end

    end

    %     if Estrella_identification == 1
    %         pos_FEDFUNDS = find(ismember(var_mnemonic_i, {'FEDFUNDS'}));
    %         pos_CPIAUCSL = find(ismember(var_mnemonic_i, {'CPIAUCSL'}));
    %         invVtheta{pos_CPIAUCSL}(pos_FEDFUNDS+1,pos_FEDFUNDS+1) = 10^9;
    %     end

    % Sigma = @inverse(L)*@makediagonal(sig)*@transpose(@inverse(L))
    if contains(prior_type,'srp','ignorecase',true)
        if hasELBdata
            % Storage Matrix
            shadowrate_all  = NaN(nloop-burnin,Nshadowrates,elbT);
            shadowrate      = NaN(Nshadowrates,elbT);
        end
    end

    % Storage Matrices
    store_Sig   = zeros(nloop-burnin,N);
    store_beta  = zeros(nloop-burnin,size(bols(:),1));
    store_L     = zeros(nloop-burnin,N*N);

    Sig = ones(N,1);
    L = eye(N);

    % Stochastic Volatility starting values
    Sigh = .05*ones(N,1);
    h0 = log(var(Y))';
    h = repmat(h0',T,1);

    if contains(prior_type,'stvol','ignorecase',true)
        store_h     = zeros(nloop-burnin,T*N);
        store_Sigh  = zeros(nloop-burnin,N);
        store_h0    = zeros(nloop-burnin,N);
    end

    if contains(prior_type,'ssp','ignorecase',true)
        store_MU = zeros(nloop-burnin,N);
    end

    %% create X block for each equation
    X_actual = X;
    % Xshadow = NaN(T,K);
    NactualBlocks = sum(actualrateBlock);
    NshadowBlocks = N - NactualBlocks;

    XX      = NaN(T,K,N);
    XX(:,:,actualrateBlock)  = repmat(X_actual, [1 1 NactualBlocks]);
    XX(:,:,~actualrateBlock) = repmat(X_actual, [1 1 NshadowBlocks]);

    YY      = NaN(T,N,N);
    YY(:,:,actualrateBlock)  = repmat(Y, [1 1 NactualBlocks]);
    YY(:,:,~actualrateBlock) = repmat(Y, [1 1 NshadowBlocks]);

    X_actual_backup = X_actual;
    Y_actual_backup = Y;

    prevdraw.Y_shadow = Y;
    prevdraw.X_shadow = X;

    beta_store_mat      = bols;

    sqrtht = exp(h/2);

    countELBaccept       = 0;
    countELBacceptBurnin = 0;
    stackAccept = NaN(nloop,1);

    for loop = 1:nloop

        if contains(prior_type,'ssp','ignorecase',true) && contains(prior_type,'srp','ignorecase',true)

            % Not demeaned data for the Y1 calculation
            dataT_shadow_not_demeaned = [Y0;prevdraw.Y_shadow];
            Y_shadow_not_demeaned = dataT_shadow_not_demeaned(p+1:size(dataT_shadow_not_demeaned,1),:);

            X_shadow_not_demeaned = zeros(T,N*p);
            tempY = dataT_shadow_not_demeaned;
            for i=1:p
                X_shadow_not_demeaned(:,(i-1)*N+1:i*N) = tempY(p-i+1:end-i,:);
            end

            dataT_actual_not_demeaned = [Y0;Y_actual_backup];
            Y_actual_not_demeaned = dataT_actual_not_demeaned(p+1:size(dataT_actual_not_demeaned,1),:);

            X_actual_not_demeaned = zeros(T,N*p);
            tempY = dataT_actual_not_demeaned;
            for i=1:p
                X_actual_not_demeaned(:,(i-1)*N+1:i*N) = tempY(p-i+1:end-i,:);
            end

            % Demeaned data for all other calculations
            dataT_shadow_demeaned = [Y0;prevdraw.Y_shadow]-repmat(MU(1:N)',size(Y0,1)+size(Y,1),1);
            Y_shadow_demeaned = dataT_shadow_demeaned(p+1:size(dataT_shadow_demeaned,1),:);

            X_shadow_demeaned = zeros(T,N*p);
            tempY = dataT_shadow_demeaned;
            for i=1:p
                X_shadow_demeaned(:,(i-1)*N+1:i*N) = tempY(p-i+1:end-i,:);
            end

            dataT_actual_demeaned = [Y0;Y_actual_backup]-repmat(MU(1:N)',size(Y0,1)+size(Y,1),1);
            Y_actual_demeaned = dataT_actual_demeaned(p+1:size(dataT_actual_demeaned,1),:);

            X_actual_demeaned = zeros(T,N*p);
            tempY = dataT_actual_demeaned;
            for i=1:p
                X_actual_demeaned(:,(i-1)*N+1:i*N) = tempY(p-i+1:end-i,:);
            end

            %% create X block for each equation
            NactualBlocks = sum(actualrateBlock);
            NshadowBlocks = N - NactualBlocks;

            XX                       = NaN(T,K,N);
            XX(:,:,actualrateBlock)  = repmat(X_actual_demeaned, [1 1 NactualBlocks]);
            XX(:,:,~actualrateBlock) = repmat(X_shadow_demeaned, [1 1 NshadowBlocks]);

            YY                       = NaN(T,N,N);
            YY(:,:,actualrateBlock)  = repmat(Y_actual_demeaned, [1 1 NactualBlocks]);
            YY(:,:,~actualrateBlock) = repmat(Y_shadow_demeaned, [1 1 NshadowBlocks]);

            prevdraw.MU         = MU(1:N)';

        end

        if contains(prior_type,'ssp','ignorecase',true) && ~contains(prior_type,'srp','ignorecase',true)
            dataT_demeaned = [Y0;Y_actual_backup]-repmat(MU(1:N)',size(Y0,1)+size(Y,1),1);
            Y_demeaned = dataT_demeaned(p+1:size(dataT_demeaned,1),:);
            X_demeaned = zeros(T,N*p);
            tempY = dataT_demeaned;
            for i=1:p
                X_demeaned(:,(i-1)*N+1:i*N) = tempY(p-i+1:end-i,:);
            end

            %% create X, Y blocks for each equation
            NactualBlocks = sum(actualrateBlock);
            NshadowBlocks = N - NactualBlocks;

            XX      = NaN(T,K,N);
            XX(:,:,actualrateBlock)  = repmat(X_demeaned, [1 1 NactualBlocks]);
            XX(:,:,~actualrateBlock) = repmat(X_demeaned, [1 1 NshadowBlocks]);

            YY      = NaN(T,N,N);
            YY(:,:,actualrateBlock)  = repmat(Y_demeaned, [1 1 NactualBlocks]);
            YY(:,:,~actualrateBlock) = repmat(Y_demeaned, [1 1 NshadowBlocks]);

            prevdraw.MU         = MU(1:N)';
        end

        if ~contains(prior_type,'ssp','ignorecase',true) && contains(prior_type,'srp','ignorecase',true)

            dataT_shadow = [Y0;prevdraw.Y_shadow];
            Y_shadow = dataT_shadow(p+1:size(dataT_shadow,1),:);

            X_shadow = zeros(T,N*p);
            tempY = dataT_shadow;
            for i=1:p
                X_shadow(:,(i-1)*N+1:i*N) = tempY(p-i+1:end-i,:);
            end
            X_shadow                        = [ones(T,1) X_shadow];

            dataT_actual = [Y0;Y_actual_backup];
            Y_actual = dataT_actual(p+1:size(dataT_actual,1),:);

            X_actual = zeros(T,N*p);
            tempY = dataT_actual;
            for i=1:p
                X_actual(:,(i-1)*N+1:i*N) = tempY(p-i+1:end-i,:);
            end
            X_actual                        = [ones(T,1) X_actual];

            %% create X block for each equation
            NactualBlocks = sum(actualrateBlock);
            NshadowBlocks = N - NactualBlocks;

            XX                       = NaN(T,K,N);
            XX(:,:,actualrateBlock)  = repmat(X_actual, [1 1 NactualBlocks]);
            XX(:,:,~actualrateBlock) = repmat(X_shadow, [1 1 NshadowBlocks]);

            YY                       = NaN(T,N,N);
            YY(:,:,actualrateBlock)  = repmat(Y_actual, [1 1 NactualBlocks]);
            YY(:,:,~actualrateBlock) = repmat(Y_shadow, [1 1 NshadowBlocks]);

        end

        stationarity = 0;
        while stationarity == 0 % truncate non-stationary draws

            for kh_i = 1:thin_gibbs

                for ii = 1:N

                    % Start Simulation
                    if ii == 1
                        Y = YY(:,:,ii);
                        zt = XX(:,:,ii);
                        km = k/N;
                        unified_simulation
                    else
                        Y = YY(:,:,ii);
                        zt = [XX(:,:,ii) -YY(:,1:ii-1,ii)];
                        km = k/N + ii-1;
                        unified_simulation
                    end
                end

                % Store results
                beta_store = [];
                a_store    = [];

                %compute sqrtht^0.5 from volatility states, needed in shadow rate step
                if contains(prior_type,'stvol','ignorecase',true)
                    sqrtht     = exp(h/2);
                else
                    h          = repmat(log(Sig)',T,1);
                    sqrtht     = exp(h/2);
                end

                for iv = 1:N
                    if iv ==1
                        km = k/N;
                        theta1 = theta{1,1};
                        beta_store =[beta_store;theta1];
                    else
                        km = k/N + iv-1;
                        theta1 = theta{iv,1};
                        beta_store =[beta_store;theta1(1:k/N)];
                        a_store =[a_store;theta1(k/N+1:end)];
                        L(iv,1:iv-1) = theta1(k/N+1:end)';
                    end
                end

                % Convert structural parameters to reduce form parameters
                beta_store_mat_i = reshape(beta_store,K,N)';
                if contains(prior_type,'ssp','ignorecase',true)
                    beta_store_mat = [L\beta_store_mat_i(:,1:N) L\beta_store_mat_i(:,N+1:2*N)]';
                else
                    %                         beta_store_mat = [L\beta_store_mat_i(:,1) L\beta_store_mat_i(:,2:N+1) L\beta_store_mat_i(:,N+2:2*N+1)]';
                    beta_store_mat = L\beta_store_mat_i(:,1);
                    stepi=0;
                    for ki = 1:p
                        beta_store_mat = [beta_store_mat L\beta_store_mat_i(:,2+stepi:N+1+stepi)];
                        stepi = stepi+N;
                    end
                    beta_store_mat = beta_store_mat';
                end

                beta_store_structural = beta_store;
                beta_store=beta_store_mat(:);

                % DL: Dirichlet-Laplace
                if contains(prior_type,'DirLap','ignorecase',true)
                    %                     invVtheta_backup=invVtheta;

                    %% sample psibeta
                    nupsibeta = varthetabeta.*(zetabeta./abs(beta_store_structural));
                    invpsibeta = random('inversegaussian',nupsibeta,1);
                    psibeta = 1./(invpsibeta  + 10^-6);
                    %% sample zetabeta
                    zetabeta = gigrnd((k)*(abeta-1),1,2*sum(abs(beta_store_structural)./varthetabeta),1);
                    %% sample varthetabeta
                    bigLbeta = zeros(k,1);
                    for tv = 1:k
                        bigLbeta(tv,:) = gigrnd(abeta-1,1,2*abs(beta_store_structural(tv,:)),1);
                    end
                    varthetabeta = bigLbeta/sum(bigLbeta);
                    if contains(prior_type,'ssp','ignorecase',true)
                        invVbeta = reshape(1./(psibeta.*(varthetabeta.^2).*zetabeta.^2 + 10^-6),np,N)';
                    else
                        invVbeta = reshape(1./(psibeta.*(varthetabeta.^2).*zetabeta.^2 + 10^-6),np+1,N)';
                    end

                    %% sample psia
                    nupsia = varthetaa.*(zetaa./abs(a_store));
                    invpsia = random('inversegaussian',nupsia,1);
                    psia = 1./(invpsia  + 10^-6);
                    %% sample zetaa
                    zetaa = gigrnd((m)*(agam-1),1,2*sum(abs(a_store)./varthetaa),1)  + 10^-6;
                    %% sample varthetaa
                    bigLa = zeros(m,1);
                    for tv = 1:m
                        bigLa(tv,:) = gigrnd(agam-1,1,2*abs(a_store(tv,:)),1);
                    end
                    varthetaa = bigLa/sum(bigLa )  + 10^-6;
                    invVa = 1./(psia.*(varthetaa.^2).*zetaa.^2 + 10^-6);

                    for iv = 1:N
                        invVtheta{iv,1} = sparse(diag([invVbeta(iv,:)';invVa(idn==iv)]));
                    end

                end

                % Steady State Prior
                if contains(prior_type,'ssp','ignorecase',true)

                    beta_store_mat_sans_constant = reshape(beta_store,K,N);

                    Y1=[];
                    if contains(prior_type,'srp','ignorecase',true)
                        for var_i = 1:N
                            %                                 Y1_i = (YY(:,:,var_i) + prevdraw.MU)-(XX(:,:,var_i)+repmat(prevdraw.MU,T,2))*beta_store_mat_sans_constant(:,var_i);
                            if actualrateBlock(var_i) == 1
                                Y1_i = Y_actual_not_demeaned(:,var_i) - X_actual_not_demeaned*beta_store_mat_sans_constant(:,var_i);
                            else
                                Y1_i = Y_shadow_not_demeaned(:,var_i)-X_shadow_not_demeaned*beta_store_mat_sans_constant(:,var_i);
                            end

                            Y1(:,var_i) = Y1_i;
                        end
                    else
                        Y1 = Y_actual_backup-X_actual_backup*beta_store_mat_sans_constant;
                    end

                    U = eye(N);
                    jj = 1;
                    for jx = 1:p
                        betai = beta_store_mat_sans_constant(jj:jj+N-1,:);
                        U=[U;betai'];
                        jj=jj+N;
                    end
                    D = [ones(T,1) -ones(T,p)];
                    if contains(prior_type,'stvol','ignorecase',true)
                        sigma_1                     = (L\eyeN)*diag(exp(rowextract(h,1)))*(L\eyeN)';
                        invSigma_1                  = sigma_1\eyeN;
                        Sum_DD_invSigma             = kron(D'*D,invSigma_1);
                        vec_Sum_invSigma_Y1_D       = vec(invSigma_1*Y1'*D);
                        for tt = 2:T
                            sigma_1                 = (L\eyeN)*diag(exp(rowextract(h,1)))*(L\eyeN)';
                            invSigma_tt             = sigma_1\eyeN;
                            %                                 invSigma_tt             = L'*diag(exp(rowextract(-h,tt)))*L;
                            Sum_DD_invSigma         = Sum_DD_invSigma + kron(D'*D,invSigma_tt);
                            vec_Sum_invSigma_Y1_D   = vec_Sum_invSigma_Y1_D + vec(invSigma_tt*Y1'*D);
                        end
                        Vstar                       = (invV0 + U'*Sum_DD_invSigma*U)\eyeN;
                        MUstar                      = Vstar*(invV0*M0 + U'*vec_Sum_invSigma_Y1_D);
                        MU                          = MUstar + (randn(1,N)*chol(Vstar + 1e-13*eye(size(Vstar))))';
                    else
                        sigma_1                     = (L\eyeN)*diag(Sig)*(L\eyeN)';
                        invSigma                    = sigma_1\eyeN;
                        Vstar                       = (invV0 + U'*kron(D'*D,invSigma)*U)\eyeN;
                        MUstar                      = Vstar*(invV0*M0+U'*vec(invSigma*Y1'*D));
                        MU                          = MUstar + (randn(1,N)*chol(Vstar + 1e-13*eye(size(Vstar))))';
                    end

                    if contains(prior_type,'srp','ignorecase',true)
                        % Bound MU to the ELB
                        ndxMUshadow=(MU(ndxSHADOWRATE)<ELBbound');
                        MU(ndxSHADOWRATE)=ELBbound'.*ndxMUshadow + MU(ndxSHADOWRATE).*(1-ndxMUshadow);
                    end

                    if SSPssvs == 1
                        prob_MU           = (normpdf(M0,MU,sqrt(tau1_MU))*(1-probability_i_MU))./((normpdf(M0,MU,sqrt(tau1_MU))*(1-probability_i_MU)) + ...
                            (normpdf(M0,MU,sqrt(tau0_MU))*probability_i_MU) );
                        gam1_MU            = prob_MU > rand( size(prob_MU,1),1);
                        DD1_MU           = (gam1_MU).*tau1_MU + (1-gam1_MU).*tau0_MU;
                        invV0 = sparse(diag(1./DD1_MU));
                    end

                    % add the simulated / implied constant to the coefficient matrix
                    F = [beta_store_mat_sans_constant' ; eye(N*(p-1),N*p)];
                    Csim = (eye(size(F,1)) -F)*repmat(MU,p,1);
                    beta_store_mat = [Csim(1:N)'; beta_store_mat_sans_constant];

                else

                    beta_store_mat = reshape(beta_store,K,N);

                end % if contains(prior_type,'ssp'

                if check_stationarity == 1

                    F_comp = [beta_store_mat(2:end,:)' ; [eye(N*(p-1)),zeros(N*(p-1),N)]];
                    if max(abs(eig(F_comp)))>1
                        stationarity = 0;
                    else
                        stationarity = 1;
                    end

                else
                    stationarity = 1;
                end

            end

        end % stationarity

        % data.*repmat(sigma_unscale,T+p,1)+repmat(mi_unscale,T+p,1)+repmat(oldMU(1:N)',T+p,1)
        if contains(prior_type,'srp','ignorecase',true)

            if hasELBdata

                if contains(prior_type,'ssp','ignorecase',true)
                    %                         Y = Y_shadow_demeaned + prevdraw.MU;
                    Y = Y_shadow_not_demeaned;
                else
                    Y = Y_shadow;
                end

                % update state space objects
                elb.Y               = Y(elb.T0+1:end,:)';

                % adjust obs for lagged actual rates
                beta_store_mat_actual                     = beta_store_mat(ndxSHADOWRATELAGS,:);
                beta_store_mat_actual(:,~actualrateBlock) = 0;
%                 Yhatactual = Xactual(elb.T0+1:end,ndxSHADOWRATELAGS) * beta_store_mat_actual;
                Yhatactual = transpose(Xactual(elbT0+1:end,ndxSHADOWRATELAGS) * beta_store_mat_actual);
                %                 elb.Y      = elb.Y - Yhatactual';

                % update VAR companion form
                beta_store_mat_shadow                                    = beta_store_mat;
                beta_store_mat_shadow(ndxSHADOWRATELAGS,actualrateBlock) = 0;
                elb.A(elb.ndxY, :)                                       = beta_store_mat_shadow';

                elb.B(elb.ndxY,:)   = L\eyeN;

                % map SV into sqrtSigma
                elb.sqrtSigma = sqrtht(elb.T0+1:end,:)';

                elbYgibbs           = elb.Y;
                elb.Y(elb.yNaN)     = NaN; % missing values

                %% PS setup
                beta_store_mat0    = transpose(beta_store_mat_shadow(1,:)) + Yhatactual;
                beta_store_mat3    = reshape(transpose(beta_store_mat_shadow(2:end,:)), N, N, p);
                invbbb  = L ./ permute(elb.sqrtSigma, [1 3 2]);
                elbY0   = reshape(elb.X0(2:end), N, p);
                if loop == 1
                    [~, CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx] = VARTVPSVprecisionsamplerNaN(beta_store_mat3,invbbb,elb.Y,elb.yNaN,elbY0,beta_store_mat0);
                end

                % b) reconstruct Y and X
                shadowYdata = Ydata;

                %
                %                 shadowrate_draws = gibbsdrawShadowrates(elb.Y, elb.X0, elb.ndxS, elb.sNaN, p, elb.A, elb.B, elb.sqrtSigma, ...
                %                     elb.bound, 1, elb.gibbsburn);
                %
                %                 shadowrate = mean(shadowrate_draws,3);

                if loop < burnin * .5
                    shadowrate = gibbsdrawShadowrates(elbYgibbs, elb.X0, Yhatactual, elb.ndxS, elb.sNaN, p, elb.A, elb.B, elb.sqrtSigma, ...
                        elb.bound, 1, elb.gibbsburn);
                else % try acceptance sampling
                    YYdraws = VARTVPSVprecisionsamplerNaN(beta_store_mat3,invbbb,elb.Y,elb.yNaN,...
                        elbY0,beta_store_mat0, ...
                        CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx, elb.Nproposals);
                    YYdraws = reshape(YYdraws, N, elbT, elb.Nproposals);

                    shadowrateProposals = YYdraws(ndxSHADOWRATE,:,:);

                    Ok = false;
                    ndxAccept = 0;
                    while ~Ok && (ndxAccept < elb.Nproposals)
                        ndxAccept = ndxAccept + 1;
                        thisProposal = shadowrateProposals(:,:,ndxAccept);
                        Ok = all(thisProposal(elb.sNaN) < ELBbound);
                    end
                    if Ok
                        shadowrate     = shadowrateProposals(:,:,ndxAccept);
                        if loop > burnin
                            countELBaccept = countELBaccept + 1;
                            stackAccept(loop-burnin) = ndxAccept;
                        else
                            countELBacceptBurnin = countELBacceptBurnin + 1;
                        end
                    else
                        shadowrate = gibbsdrawShadowrates(elbYgibbs, elb.X0, Yhatactual, elb.ndxS, elb.sNaN, p, elb.A, elb.B, elb.sqrtSigma, ...
                            elb.bound, 1, elb.gibbsburn);
                        % fprintf('%d, none accepted\n', m);
                    end
                end

%                 shadowYdata(p+elb.T0+1:end,ndxSHADOWRATE) = shadowrate';
                  shadowYdata(p+elb.T0+1:end,ndxSHADOWRATE) = shadowYdata(p+elb.T0+1:end,ndxSHADOWRATE).*(1-ELBdummy(p+elb.T0+1:end,:)) + ...
                                                              shadowrate'.*ELBdummy(p+elb.T0+1:end,:);

                % matrix X
                lags=zeros(Nobs,N*p);
                for l=1:p
                    lags(p+1:Nobs,(N*(l-1)+1):N*l) = shadowYdata(p+1-l:Nobs-l,1:N);
                end

                if contains(prior_type,'ssp','ignorecase',true)
                    X = [lags(p+1:Nobs,:)];
                else
                    X = [ones(Nobs-p,1) lags(p+1:Nobs,:)];
                end
                Y = shadowYdata(p+1:end,:);

                if NORMALIZE_DATA==1
                    Y_original_shadow = Y.*sigma_unscale+mi_unscale;
                    shadowrate = Y_original_shadow(elbT0+1:end,ndxSHADOWRATE);
                else
                    shadowrate = Y(elbT0+1:end,ndxSHADOWRATE);
                end

                % generate state vector for forecast jumpoff
                Xjumpoff     = zeros(Nstates,1);
                Xjumpoff(1) = 1;
                for l=1:p
                    Xjumpoff(1+(l-1)*N+(1:N)) = Y(end-(l-1),1:N);
                end
                % add lagged actual rates to jumpoff
                for l=1:p
                    Xjumpoff(size(bols,1)+(l-1)*Nyields+(1:Nyields)) = data(Nobs-(l-1),ndxYIELDS);
                end

                prevdraw.Y_shadow = Y;
                prevdraw.X_shadow = X;

            end % hasELBdata

        end


        %% Storage
        if loop>burnin

            thisdraw = loop-burnin;

            if contains(prior_type,'stvol','ignorecase',true)
                store_h(thisdraw,:) = reshape(h',1,T*N);
                store_Sigh(thisdraw,:) = Sigh';
                store_h0(thisdraw,:) = h0';
            end

            if contains(prior_type,'ssp','ignorecase',true)
                if NORMALIZE_DATA==1
                    store_MU(thisdraw,:) = (MU(1:N,:).*sigma_unscale'+mi_unscale')';
                    %                         store_MU(thisdraw,:) = MU(1:N,:);
                else
                    store_MU(thisdraw,:) = MU(1:N,:);
                end
            end

            store_beta(thisdraw,:) = beta_store_mat(:)';
            store_L(thisdraw,:) = L(:);

            if contains(prior_type,'srp','ignorecase',true)
                if hasELBdata
                    shadowrate_all(thisdraw,:,:)  = shadowrate';
                end
            end


            % Out of Sample forecasting
            if draw_predictive_densities == 1

                % draw and scale SV shocks
                logSV0          = h(end,:)'; % Note: Vol_states record logs of *variances*

                % IRF
                stnd_norm_draws_irf = randn(N, irf_forecast_horizon * Ndraws);
                logSVshocks_irf     = sqrt(Sigh).*stnd_norm_draws_irf;
                logSVshocks_irf     = reshape(logSVshocks_irf, N, irf_forecast_horizon, Ndraws);

                logSV_irf           = logSV0 + cumsum(logSVshocks_irf,2);
                fcstSVdraws_irf     = exp(logSV_irf * 0.5);

                % draw random numbers
                nushocks_irf        = fcstSVdraws_irf .* randn(N, irf_forecast_horizon, Ndraws);

                for nn = 1 : Ndraws

                    if contains(prior_type,'srp','ignorecase',true)
                        % update VAR companion form
                        PAIactual                     = beta_store_mat(ndxYIELDLAGS,:);
                        PAIactual(:,~actualrateBlock) = 0;

                        PAIshadow = beta_store_mat;
                        PAIshadow(ndxYIELDLAGS,actualrateBlock) = 0;

                        fcstA(ndxfcstY, 1:size(PAIshadow',2))         = PAIshadow';
                        fcstA(ndxfcstY, size(PAIshadow',2)+1:Nstates) = PAIactual';

                        % forecast simulation
                        fcstX0 = Xjumpoff;

                    else
                        % update VAR companion form
                        fcstA(ndxfcstY, 1:(N*p+1)) = beta_store_mat';

                        % forecast simulation
                        Xjumpoff_backup=Xjumpoff;
                        Xjumpoff = Xjumpoff_backup(1:size(bols,1));
                        fcstX0 = Xjumpoff;
                    end

                    %% GIRFs

                    if draw_GIRFs == 1

                        GIRFs

                    end


                end % forecast Ndraws

            end % if predictive densities

        end % if loop>burnin

    end % for loop

    if draw_GIRFs == 1

        %% censor IRFsims if necessary
        % -----------------------------------------------------------------------------------------
        fcstYdraws_backup = fcstYdraws_irf;

        if contains(prior_type,'srp','ignorecase',true)
            % no shock
            yieldDraws1                         = fcstYdraws_irf(ndxYIELDS,:,:,:);
            ndx                                 = yieldDraws1 < ELBboundyields';
            yieldDraws1                         = ndx.*ELBboundyields'+(1-ndx).*yieldDraws1;
            fcstYdraws_irf(ndxYIELDS,:,:,:)     = yieldDraws1;
            %                 fcstYdraws_irf(cumcode,:,:)         = cumsum(fcstYdraws_irf(cumcode,:,:), 2);

            % plus shock
            yieldDraws2                         = fcstYdraws1plus(ndxYIELDS,:,:,:);
            ndx                                 = yieldDraws2 < ELBboundyields';
            yieldDraws2                         = ndx.*ELBboundyields'+(1-ndx).*yieldDraws2;
            fcstYdraws1plus(ndxYIELDS,:,:,:)    = yieldDraws2;
            %                 fcstYdraws1plus(cumcode,:,:)        = cumsum(fcstYdraws1plus(cumcode,:,:), 2);

            % minus shock
            yieldDraws3                         = fcstYdraws1minus(ndxYIELDS,:,:,:);
            ndx                                 = yieldDraws3 < ELBboundyields';
            yieldDraws3                         = ndx.*ELBboundyields'+(1-ndx).*yieldDraws3;
            fcstYdraws1minus(ndxYIELDS,:,:,:)   = yieldDraws3;
            %                 fcstYdraws1minus(cumcode,:,:)       = cumsum(fcstYdraws1minus(cumcode,:,:), 2);

        end

        % first, integrate only per MCMC node
        fcstYhat       = mean(fcstYdraws_irf, 3);
        fcstYhat1plus  = mean(fcstYdraws1plus, 3);
        fcstYhat1minus = mean(fcstYdraws1minus, 3);

        %% IRF
        IRFdraws1plus   = fcstYhat1plus  - fcstYhat;
        IRFdraws1minus  = fcstYhat1minus - fcstYhat;

        % integrate over MCMC nodes
        fcstYhat1plus  = mean(fcstYhat1plus, 4);
        fcstYhat1minus = mean(fcstYhat1minus, 4);

        IRF1plus        = mean(IRFdraws1plus, 4);
        IRF1plusTails   = prctile(IRFdraws1plus, prc70, 4);
        IRF1plusTails(:,:,:,3)   = median(IRFdraws1plus, 4);

        IRF1minus       = mean(IRFdraws1minus, 4);
        IRF1minusTails  = prctile(IRFdraws1minus, prc70, 4);
        IRF1minusTails(:,:,:,3)   = median(IRFdraws1minus, 4);
    end

    % collect draws
    if NORMALIZE_DATA == 1
        fcstYdraws = fcstYdraws_backup.*repmat(sigma_unscale,irf_forecast_horizon,1)'+repmat(mi_unscale,irf_forecast_horizon,1)';
    else
        fcstYdraws = fcstYdraws_backup;
    end

    % censor shadow rates only for SR models
    if contains(prior_type,'srp','ignorecase',true)
        yieldDraws                 = fcstYdraws(ndxYIELDS,:,:,:);
        ndx                        = yieldDraws < FLOOR_ELB;
        yieldDraws(ndx)            = FLOOR_ELB;
        fcstYdraws(ndxYIELDS,:,:,:)  = yieldDraws;
    end

    fcstYhat       = mean(fcstYdraws, 3);
    fcstYhat_irf_prctile = reshape(fcstYhat,N,irf_forecast_horizon,nloop-burnin);

    if contains(prior_type,'srp','ignorecase',true)
        if hasELBdata
            shadowrates     = permute(shadowrate_all, [3 2 1]);
            shadowrateMid   = median(shadowrates,3);
            shadowrateTails = prctile(shadowrates,10,3);
            shadowrateTails(:,:,2) = mean(shadowrates,3);
            shadowrateTails(:,:,3) = prctile(shadowrates,90,3);
            reshape(shadowrateTails(:,1,:),size(shadowrateTails(:,1,:),1),3);
        end
    end

    fcstYhat2       = mean(fcstYhat, 4);
    IRFPlus_tails_all_samples{vint_i,1}     = IRF1plusTails;
    IRFMinus_tails_all_samples{vint_i,1}    = IRF1minusTails;
    yforecast_all_samples{vint_i,1}         = fcstYhat2(:,forecast_horizon_set);
    fcstYdraws_all_samples{vint_i,1}        = fcstYhat_irf_prctile(:,forecast_horizon_set,:);

    if contains(prior_type,'srp','ignorecase',true)
        if hasELBdata
            shadowrateTails_all_samples{vint_i,1}  = shadowrateTails;
        end
    end

    if contains(prior_type,'ssp','ignorecase',true)
        MU_all_samples{vint_i,1}  = prctile(store_MU,[10 50 90],1)';
    end

    if contains(prior_type,'stvol','ignorecase',true)
        post_h_all_samples{vint_i,1} = reshape(mean(store_h,1),N,T)';
    end

end
