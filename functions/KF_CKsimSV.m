function [aSim]=KF_CKsimSV(yStar,Zt,R,Q,mu0,P0)
%% Purpose
%  This function runs the Kalman Filter and Backward Simulation recursions
%  from the Carter Kohn algorithm.#
% mu0 = startMean
% P0 = startVar
%% Output
%  - aSim - 
%
% 
% yStar=y_KF;
% Zt=Z_KF;
% R=R_KF;
% Q=eye(N)*0.0001;
% mu0=M0;
% P0=0.0001*eye(N);

%% Initialize
[N, T]=size(yStar);

mu_Prediction=zeros(N,T);
P_Prediction=zeros(N,N,T);
mu_Update=zeros(N,T);
P_Update=zeros(N,N,T);
%% Kalman Filter Recursion Prediction
% Following Kim Nelson (1998):
for t=1:T    
    % Prediction
    if t==1
        mu_Prediction(:,t) = mu0; %at+1|t again using the fact that T is identity
        P_Prediction(:,:,t)=P0; %variance of at+1|t using the fact that R is identity and phi = Q
    else
        mu_Prediction(:,t) = mu_Update(:,t-1); %at+1|t again using the fact that T is identity
        P_Prediction(:,:,t)=P_Update(:,:,t-1)+Q; %variance of at+1|t using the fact that R is identity and phi = Q
    end
    etat=yStar(:,t)-Zt*mu_Prediction(:,t);  % eta_{t|t-1} = y_{t} - y_{t-1} = y_{t} - x_{t}*beta_{t|t-1}
    Ft=Zt*P_Prediction(:,:,t)*Zt'+R(:,:,t); % f_{t|t-1} = x_{t}*P_{t|t-1}*x_{t}' + R
    Kt=(P_Prediction(:,:,t)*Zt')*inv(Ft); %Kalman Gain

    % Updating
    mu_Update(:,t)=mu_Prediction(:,t)+Kt*etat; % since Tt=identity                % Kim Nelson (1998): P_{t|t} = P_{t|t-1} - K_{t}*x_{t}*P_{t|t-1}
    P_Update(:,:,t)=P_Prediction(:,:,t)-Kt*Zt*P_Prediction(:,:,t);  %Pt|t, variance of the state at given Yt 
end
% last round
% etaT=yStar(:,T)-mu_Prediction(:,T);
% 
% FT=Zt*P_Prediction(:,:,T)*Zt'+R(:,:,T);
% 
% % Finv=1/FT;
% % Finv = inv(Ft);
% KT=P_Prediction(:,:,T)*Zt'*inv(FT);
% mu_Update(:,T)=mu_Prediction(:,T)+KT*etaT; % since T_T=1
% P_Update(:,:,T)=P_Prediction(:,:,T)-Kt*Zt*P_Prediction(:,:,T); % P_Prediction(:,:,T)*(1-KT);
%% Continue with the Carter Kohn (simulation) recursion %%
drawsAlpha=zeros(N,T);
% startround (t=T)
meanAlpha=mu_Update(:,T);
varAlpha=P_Update(:,:,T);
drawsAlpha(:,T)=meanAlpha+chol(varAlpha)'*randn(N,1); %draw the last observation from N(meanAlpha(T), varAlpha(T))
% iteratively obtain the draws
for t=T-1:-1:1
    % meanAlpha=mu_Update(:,t)+(P_Update(:,:,t)*inv(P_Prediction(:,:,t+1)))*(drawsAlpha(:,t+1)-mu_Update(:,t)); %compute the conditional mean of at|Yn
    % varAlpha=P_Update(:,:,t)*(eye(N)-(P_Update(:,:,t)*inv(P_Prediction(:,:,t+1)))); %compute the conditional variance of at|Yn
    % drawsAlpha(:,t)=meanAlpha+chol(varAlpha)'*randn(N,1); %compute a draw from the normal distribution for at|Yn from N(meanAlpha,varianceAlpha)

    Jt = P_Update(:,:,t) * inv(P_Prediction(:,:,t+1));
    meanAlpha = mu_Update(:,t) + Jt * (drawsAlpha(:,t+1)-mu_Update(:,t));
    varAlpha = P_Update(:,:,t) - Jt * P_Prediction(:,:,t+1) * Jt';
    drawsAlpha(:,t) = meanAlpha + chol(varAlpha)'*randn(N,1); 

end 
    
%% save the draws
aSim=drawsAlpha;
end
