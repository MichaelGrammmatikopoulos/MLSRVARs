
%% If the model has stochastic volatility sample from here
if contains(prior_type,'stvol','ignorecase',true) %#ok<STRIFCND>
    
    %% sample theta
    expinvh     = sparse(1:T,1:T,exp(-h(:,ii)),T,T);
    Ktheta      = invVtheta{ii,1} + zt'*expinvh*zt + 1e-6*eye(size(invVtheta{ii,1})); 
    thetahat    = Ktheta\(invVtheta{ii,1}*theta0{ii,1}+zt'*expinvh*Y(:,ii));
    theta{ii,1} = thetahat + chol(Ktheta,'lower')'\randn(size(Ktheta,1),1);

    % Ktheta      = nearestSPD(invVtheta{ii,1} + zt'*expinvh*zt + 10^-6*eye(size(invVtheta{ii,1}))); % 
    % thetahat    = pinv(Ktheta)*(invVtheta{ii,1}*theta0{ii,1}+zt'*expinvh*Y(:,ii));
    % theta{ii,1} = thetahat + pinv(chol(Ktheta,'lower'))'*randn(size(Ktheta,1),1);

    %% sample h
    e       = Y(:,ii) - zt*theta{ii,1};
    Ystar   = log(e.^2 + .0001);
    h(:,ii) = SVRW(Ystar,h(:,ii),Sigh(ii),h0(ii));

    %% sample h0
    Kh0    = (1./Sigh(ii) + 1./Vh);
    h0hat  = Kh0\(ah./Vh + h(1,ii)'./Sigh(ii));
    h0(ii) = h0hat + chol(Kh0,'lower')'\randn(1,1);

    %% sample Sigh
    eh                  = h(:,ii)- [h0(ii);h(1:T-1,ii)];
    Sigh(ii)            = 1./gamrnd(nu0h+T/2, 1./(S0h + sum(eh.^2)'/2));
    % if Sigh(ii)<10^-6
    %     Sigh(ii)
    %     Sigh(ii) = 10^-6;
    % end

else % For the homoskedasticity assumption sample from here
    
    %% sample theta
    Ktheta      = invVtheta{ii,1} + zt'*zt./Sig(ii) + 1e-6*eye(size(invVtheta{ii,1}));    
    thetahat    = Ktheta\(invVtheta{ii,1}*theta0{ii,1}+zt'*Y(:,ii)./Sig(ii));
    theta{ii,1} = thetahat + chol(Ktheta,'lower')'\randn(size(Ktheta,1),1);
    
    % Ktheta      = nearestSPD(invVtheta{ii,1} + zt'*zt./Sig(ii)+ 10^-6*eye(size(invVtheta{ii,1}))); %+ 10^-6*eye(size(invVtheta{ii,1}))); %+ 1e-6*eye(size(invVtheta{ii,1}));    
    % thetahat    = pinv(Ktheta)*(invVtheta{ii,1}*theta0{ii,1}+zt'*Y(:,ii)./Sig(ii));
    % theta{ii,1} = thetahat + pinv(chol(Ktheta,'lower'))'*randn(size(Ktheta,1),1);

    %% sample Sig
    e                = Y(:,ii) - zt*theta{ii,1};
    Sig(ii)          = 1./gamrnd(nu0+T/2, 1./(S0 + sum(e.^2)/2));
    % if Sig(ii)<10^-6
    %     Sigh(ii)
    %     Sig(ii) = 10^-6;
    % end
   
end

%% Machine Learning Priors for the coefficients

if contains(prior_type,'blasso','ignorecase',true) %#ok<STRIFCND>

    if contains(prior_type,'blasso_G_','ignorecase',true)  % Global shrinkage
        
        % Lambda
        lambda{ii,1}              = gamrnd(a0_global+km,1./(0.5*sum(tau{ii,1})+b0_global));

        % tau
        Stau                      = sqrt(lambda{ii,1}./(theta{ii,1}.^2));
        invtau                    = random('inversegaussian',Stau,lambda{ii,1});
        
        invtau (isnan( invtau  )) = 10^-6;
        invtau(isnan(invtau)|(invtau<=0)) = 10^-6; 
        invtau(sign(invtau).*isinf(invtau)<0) = 10^-6;
        invtau(sign(invtau).*isinf(invtau)>0) = 10^6;
        
        tau{ii,1}                 = 1./invtau;
        invVtheta{ii,1}           = (diag(invtau));

    end

    if contains(prior_type,'blasso_A_','ignorecase',true)  % Adaptive shrinkage
%         lambda{ii,1}              = gamrnd(a0_adaptive+1,1./(0.5*(tau{ii,1})+b0_adaptive));
    
        % Lambda
        lambda_c{ii,1}              = gamrnd(a0_c+1, 1./(0.5*(tau_c{ii})+b0_c));
        lambda_L{ii,1}              = gamrnd(a0_L+1, 1./(0.5*(tau_L{ii})+b0_L));
        lambda{ii,1}                = [lambda_c{ii,1};lambda_L{ii,1}];
% lambda{ii,1};
        % tau
        Stau_c                          = sqrt(lambda_c{ii,1}./((theta{ii,1}(1:no_of_regression_coef,:)).^2));
        Stau_L                          = sqrt(lambda_L{ii,1}./((theta{ii,1}(no_of_regression_coef+1:end,:)).^2));
        Stau                            = [Stau_c;Stau_L];
        Stau(sign(Stau).*isinf(Stau)<0) = 10^-6;
        Stau(sign(Stau).*isinf(Stau)>0) = 10^6;

        invtau                          = random('inversegaussian',Stau,lambda{ii,1});   
        invtau(isnan(invtau)|(invtau<=0)) = 10^-6; 
        invtau(sign(invtau).*isinf(invtau)<0) = 10^-6;
        invtau(sign(invtau).*isinf(invtau)>0) = 10^6;

        tau{ii}                     = 1./invtau;

        tau_c{ii}                   = tau{ii}(1:no_of_regression_coef,:);
        tau_L{ii}                   = tau{ii}(no_of_regression_coef+1:end,:);

        invVtheta{ii,1}             = sparse(diag(invtau));

    end


end

if contains(prior_type,'ssvs','ignorecase',true)

    % SVSS
    theta1_i        = theta{ii,1};
    theta0_i        = theta0{ii,1};

    if contains(prior_type,'ssvs_Min','ignorecase',true) % SVSS - Minnesota
        tau1_i = (tau1_min./(invVtheta_min{ii}));
        tau0_i = (tau0_min./(invVtheta_min{ii}));
        prob1       = (normpdf(theta0_i,theta1_i,sqrt(tau1_i))*(1-probability_i))./((normpdf(theta0_i,theta1_i,sqrt(tau1_i))*(1-probability_i)) + ...
            (normpdf(theta0_i,theta1_i,sqrt(tau0_i))*probability_i) );
        gam1            = prob1 > rand(size(prob1,1),1);
        DD1             = (gam1).*tau1_i + (1-gam1).*tau0_i;
        invVtheta{ii,1} = (diag(1./DD1));

    else
        tau0_i = tau0_L*ones(km,1);
        tau0_i(1:K) = tau0_c;
        
        tau1_i = tau1_L*ones(km,1);
        tau1_i(1:K) = tau1_c;

        prob1       = (normpdf(theta0_i,theta1_i,sqrt(tau1_i))*(1-probability_i))./((normpdf(theta0_i,theta1_i,sqrt(tau1_i))*(1-probability_i)) + ...
            (normpdf(theta0_i,theta1_i,sqrt(tau0_i))*probability_i) );
        gam1            = prob1 > rand(size(prob1,1),1);
        DD1             = (gam1).*tau1_i + (1-gam1).*tau0_i;
        invVtheta{ii,1} = (diag(1./DD1));
    end

end

        % invtau (isnan( invtau  )) = 10^-4;
        % if sum(invtau<1e-6)>0 
        %     1./invtau'
        % end
        % if sum(invtau<0)>0
        %     1./invtau'
        %          error("impossible negative number")
        % end
        % 
% If the model has SSVS prior sample from here
% if contains(prior_type,'ssvs','ignorecase',true)
%
%     % SVSS - Minnesota
%     if contains(prior_type,'ssvs_Min','ignorecase',true)
%         % SVSS
%         theta1_i        = theta{ii,1};
%         theta0_i        = theta0{ii,1};
%
%         tau1_coef_i = (tau1_coef./(invVtheta_min{ii}(1:no_of_regression_coef,:)));
%         tau0_coef_i = (tau0_coef./(invVtheta_min{ii}(1:no_of_regression_coef,:)));
%
%         tau1_L_i = (tau1_L./(invVtheta_min{ii}(no_of_regression_coef+1:end,:)));
%         tau0_L_i = (tau0_L./(invVtheta_min{ii}(no_of_regression_coef+1:end,:)));
%
%         % Modify such that covariances can get different shrinkage
%         prob1_coef           = (normpdf(theta0_i(1:no_of_regression_coef,:),theta1_i(1:no_of_regression_coef,:),sqrt(tau1_coef_i))*(1-probability_i_coef))./...
%             ((normpdf(theta0_i(1:no_of_regression_coef,:),theta1_i(1:no_of_regression_coef,:),sqrt(tau1_coef_i))*(1-probability_i_coef)) + ...
%             (normpdf(theta0_i(1:no_of_regression_coef,:),theta1_i(1:no_of_regression_coef,:),sqrt(tau0_coef_i))*probability_i_coef) );
%
%         prob1_L           = (normpdf(theta0_i(no_of_regression_coef+1:end,:),theta1_i(no_of_regression_coef+1:end,:),sqrt(tau1_L_i))*(1-probability_i_L))./...
%             ((normpdf(theta0_i(no_of_regression_coef+1:end,:),theta1_i(no_of_regression_coef+1:end,:),sqrt(tau1_L_i))*(1-probability_i_L)) + ...
%             (normpdf(theta0_i(no_of_regression_coef+1:end,:),theta1_i(no_of_regression_coef+1:end,:),sqrt(tau0_L_i))*probability_i_L) );
%
%         gam1_coef         = prob1_coef > rand(size(prob1_coef,1),1);
%         gam1_L            = prob1_L    > rand(size(prob1_L,1),1);
%
%         DD1_coef = gam1_coef.*tau1_coef_i + (1-gam1_coef).*tau0_coef_i;
%         DD1_L    = gam1_L.*tau1_L_i + (1-gam1_L).*tau0_L_i;
%
%         DD1 = [DD1_coef;DD1_L];
%
%         invVtheta{ii,1} = sparse(diag(1./DD1));
%
%     else
%         % SVSS
%         theta1_i        = theta{ii,1};
%         theta0_i        = theta0{ii,1};
%
%         prob1           = (normpdf(theta0_i,theta1_i,sqrt(tau1))*(1-probability_i))./((normpdf(theta0_i,theta1_i,sqrt(tau1))*(1-probability_i)) + ...
%             (normpdf(theta0_i,theta1_i,sqrt(tau0))*probability_i) );
%         gam1            = prob1 > rand(size(prob1,1),1);
%         DD1             = (gam1).*tau1 + (1-gam1).*tau0;
%         invVtheta{ii,1} = sparse(diag(1./DD1));
%
%         % SVSS
% %         theta1_i        = theta{ii,1};
% %         theta0_i        = theta0{ii,1};
% %
% %         % Modify such that covariances can get different shrinkage
% %         prob1_coef           = (normpdf(theta0_i(1:no_of_regression_coef,:),theta1_i(1:no_of_regression_coef,:),sqrt(tau1_coef_pure))*(1-probability_i_coef))./...
% %             ((normpdf(theta0_i(1:no_of_regression_coef,:),theta1_i(1:no_of_regression_coef,:),sqrt(tau1_coef_pure))*(1-probability_i_coef)) + ...
% %             (normpdf(theta0_i(1:no_of_regression_coef,:),theta1_i(1:no_of_regression_coef,:),sqrt(tau0_coef_pure))*probability_i_coef) );
% %
% %         prob1_L           = (normpdf(theta0_i(no_of_regression_coef+1:end,:),theta1_i(no_of_regression_coef+1:end,:),sqrt(tau1_L_pure))*(1-probability_i_L))./...
% %             ((normpdf(theta0_i(no_of_regression_coef+1:end,:),theta1_i(no_of_regression_coef+1:end,:),sqrt(tau1_L_pure))*(1-probability_i_L)) + ...
% %             (normpdf(theta0_i(no_of_regression_coef+1:end,:),theta1_i(no_of_regression_coef+1:end,:),sqrt(tau0_L_pure))*probability_i_L) );
% %
% %         gam1_coef         = prob1_coef > rand(size(prob1_coef,1),1);
% %         gam1_L            = prob1_L    > rand(size(prob1_L,1),1);
% %
% %         DD1_coef = gam1_coef*tau1_coef_pure + (1-gam1_coef).*tau0_coef_pure;
% %         DD1_L    = gam1_L*tau1_L_pure + (1-gam1_L).*tau0_L_pure;
% %
% %         DD1 = [DD1_coef;DD1_L];
% %
% %         invVtheta{ii,1} = sparse(diag(1./DD1));
%
%     end
%
%     id_ssvs=1;
%
% end


% if contains(prior_type,'blasso','ignorecase',true) %#ok<STRIFCND>
%
%     % %% Lambda
%     % lambda{ii,1}              = gamrnd(a0+1,1./(0.5*(tau{ii,1})+b0));
%     % % lambda{ii,1}              = gamrnd(a0+km,1./(0.5*sum(tau{ii,1})+b0));
%     %
%     % %% tau
%     % Stau                      = sqrt(lambda{ii,1}./(theta{ii,1}.^2));
%     % invtau                    = random('inversegaussian',Stau,lambda{ii,1});
%     % invtau (isnan( invtau  )) = 10^-6;
%     % tau{ii,1}                 = 1./invtau;
%     % invVtheta{ii,1}           = sparse(diag(invtau));
%     % id_blasso=1;
%
%     % Modify such that covariances can get different shrinkage
%
%     if contains(prior_type,'blasso_A_','ignorecase',true)
%         %% Adaptive shrinkage
%
%         % Lambda
%         lambda_coef{ii,1}           = gamrnd(a0_coef+1,1./(0.5*(tau_coef{ii,1})+b0_coef));
%         lambda_L{ii,1}              = gamrnd(a0_L+1,1./(0.5*(tau_L{ii,1})+b0_L));
%         lambda{ii,1}                = [lambda_coef{ii,1};lambda_L{ii,1}];
%
%         % tau
%         Stau_coef                   = sqrt(lambda_coef{ii,1}./(theta{ii,1}(1:no_of_regression_coef,:).^2));
%         Stau_L                      = sqrt(lambda_L{ii,1}./(theta{ii,1}(no_of_regression_coef+1:end,:).^2));
%         Stau                        = [Stau_coef;Stau_L];
%
%         invtau                    = random('inversegaussian',Stau,lambda{ii,1});
%         invtau (isnan( invtau  )) = 10^-6;
%         tau{ii,1}                 = 1./invtau;
%
%         tau_coef{ii,1}            = tau{ii,1}(1:no_of_regression_coef,:);
%         tau_L{ii,1}               = tau{ii,1}(no_of_regression_coef+1:end,:);
%
%         invVtheta{ii,1}           = sparse(diag(invtau));
%
%     else
%         %% Global shrinkage
%
%         % Lambda
%         lambda{ii,1}              = gamrnd(a0_global+km,1./(0.5*sum(tau{ii,1})+b0_global));
%
%         % tau
%         Stau                      = sqrt(lambda{ii,1}./(theta{ii,1}.^2));
%         invtau                    = random('inversegaussian',Stau,lambda{ii,1});
%         invtau (isnan( invtau  )) = 10^-6;
%         tau{ii,1}                 = 1./invtau;
%         invVtheta{ii,1}           = sparse(diag(invtau));
%     end
%
%     id_blasso=1;
%
% end
