function shadowrateDraws = gibbsdrawShadowrates(Y, STATE0, YHAT0, ndxS, ...
    sNaN, p, C, Psi, SVol, elbBound, Ndraws, burnin)
% GIBBSDRAWSHADOWRATES ...
%
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 05-Mar-2021 14:01:43 $
% $Revision : 1.00 $
% DEVELOPED : 9.9.0.1592791 (R2020b) Update 5
% FILENAME  : gibbsdrawShadowrates.m



%% process arguments


if nargin < 11 || isempty(Ndraws)
    Ndraws = 1;
end
if nargin < 12 || isempty(burnin)
    burnin = 0;
end
% if nargin < 13
%     rndStream = getDefaultStream;
% end
if nargin < 13
    showProgress = false;
end

[Ny,T] = size(Y);


S      = Y(ndxS,:);
Ns     = sum(ndxS);

shadowrateDraws = NaN(Ns,T,Ndraws);

%% prepare state space

Nstate  = Ny * p;     % without intercept
NNstate = Nstate + 1; % with    intercept

Nx            = Ny - Ns;
ndxX          = ~ndxS;
H             = zeros(Ns + Nx,NNstate);
H(:,1+(1:Ny)) = eye(Ny);

Nw  = size(Psi,2);
if Nw ~= Ny
    error('dimension mismatch');
end

PSIt = zeros(Ny, Nw, T);
psi  = Psi(1+(1:Ny),:); % w/o constant
for t = 1 : T
    PSIt(:,:,t) = psi * diag(SVol(:,t));
end

%% prepare smoothing weights
J                  = NaN(Ns, Nstate + Nx, T);
sqrtOmegaPosterior = NaN(Ns,Ns,T);

cc             = C(2:end,2:end); % drop constant
Cpowerp        = NaN(Nstate, Nstate, p+1);
Cpowerp(:,:,1) = eye(Nstate);
for k = 1 : p
    Cpowerp(:,:,k+1) = cc * Cpowerp(:,:,k);
end
Cpowerp = Cpowerp(:,1:Ny,:);


% smoothing weights
for t = 1 : T-p
    
    if any(sNaN(:,t))
        % prepare M
        M                                    = zeros(Nstate + Nw, Nw * (p + 1));
        for j = 0 : p
            M(1 : Nstate, Nw * j + (1 : Nw)) = Cpowerp(:,:,p+1-j) * PSIt(:,:,t+j);
        end
        M(Nstate + (1 : Nx), 1 : Nw)         = PSIt(ndxX,:,t);
        M(Nstate + Nx + 1 : end, 1 : Nw)     = PSIt(ndxS,:,t);
        
        
        % perform qr
        [~,R] = qr(M');
        R     = R';
        
        sqrtSigma                  = R(1:Nstate+Nx,1:Nstate+Nx);
        J(:,:,t)                   = R(Nstate+Nx+(1:Ns),1:Nstate+Nx) / sqrtSigma;
        sqrtOmegaPosterior(:,:,t)  = R(Nstate+Nx+(1:Ns), Nstate+Nx+(1:Ns));
    end
    
end

if isempty(t) % if T-p<1
    t = 0;
end

while t < T
    
    t = t + 1;
    
    if any(sNaN(:,t))
        
        k = T - t;
        Nsignal = Ny * k + Nx;
        
        thisM                               = zeros(Nsignal + Ns, Nsignal + Ns);
        for j = 0 : k
            thisM(1 : Nsignal, Nw * j + (1 : Nw)) = Cpowerp(1:Nsignal,:,k+1-j) * PSIt(:,:,t+j);
        end
        thisM(k * Ny + (1 : Nx), 1 : Nw)    = PSIt(ndxX,:,t);
        thisM(Nsignal + 1 : end, 1 : Nw)    = PSIt(ndxS,:,t);
        
        [~,R] = qr(thisM');
        R     = R';
        
        sqrtSigma                      = R(1 : Nsignal, 1 : Nsignal);
        J(:,:,t)                       = 0; % whacks out dummy values for future states used below
        J(:, Ny * (p - k) + 1 : end,t) = R(Nsignal + (1:Ns), 1 : Ny * k + Nx) / sqrtSigma; % "live values" at bottom of state vector
        sqrtOmegaPosterior(:,:,t)      = R(Nsignal + (1:Ns), Nsignal + (1:Ns));
        
        
    end
end

%% weights for truncMV in case of Ns > 1
if Ns > 1
    sqrtOmega1 = NaN(Ns, T);
    beta1      = NaN(Ns, Ns-1, T);
    for t = 1 : T
        if any(sNaN(:,t))
            vcv = sqrtOmegaPosterior(:,:,t) * sqrtOmegaPosterior(:,:,t)';
            for s = 1 : Ns
                ndxOther = (1:Ns ~= s);
                beta1(s,:,t)    = vcv(s,ndxOther) / vcv(ndxOther,ndxOther);
                sqrtOmega1(s,t) = sqrt(vcv(s,s) - vcv(s,ndxOther) / vcv(ndxOther,ndxOther) * vcv(ndxOther,s));
            end
        end
    end
    
    % checkdiff(sqrtOmega1(end,:), abs(sqrtOmegaPosterior(end,end,:)));
end

%% prepare further state-space objects

Cex1 = C(2:end,2:end);
Hex1 = eye(Ny,Nstate);


HC    = Hex1 * Cex1;
CCpp1 = Cex1^(p+1); % including constant

%% compute deterministic state
Y0 = zeros(Ny,T);
for t = 1 : T
    Y0(:,t) = H * STATE0;
    if ~isempty(YHAT0)
        Y0(:,t) = Y0(:,t) + YHAT0(:,t);
    end
    STATE0  = C * STATE0;
end
Ytilde = Y - Y0;

%% Gibbs draws
totalNdraws = burnin + Ndraws;

if isempty(elbBound)
    zdraws = randn(Ns, T, totalNdraws);
else
    udraws = rand(Ns, T, totalNdraws);
end



if showProgress
    progressbar(0)
end
for n = 1 : totalNdraws
    
    % init t = 1
    STATElag        = zeros(Nstate,1); 

    % prepare construction of STATE vector
    YY      = cat(2, Ytilde, zeros(Ny, p)); % padded (dummy) values to compute STATEfuture below
    
    for t = 1 : T
        
        if any(sNaN(:,t))
            
            Yhat            = HC * STATElag;
            Xresid          = Ytilde(ndxX,t) - Yhat(ndxX);
            
            STATEfuture     = YY(:,t+p:-1:t+1); % using zeros as dummy values for t + k > T, whacked out by J=0
            STATEfuturehat  = CCpp1 * STATElag;
            STATEtilde      = STATEfuture(:) - STATEfuturehat;

            Shat            = Yhat(ndxS) + Y0(ndxS,t);
            Sposterior      = Shat + J(:,:,t) * [STATEtilde; Xresid]; % note: could adapt to t > T - p

            if isempty(elbBound)
                S(:,t) = Sposterior + sqrtOmegaPosterior(:,:,t) * zdraws(:,t,n);
            else
                if Ns == 1
                    
                    S(:,t) = drawTruncNormal(Sposterior, sqrtOmegaPosterior(:,:,t), elbBound, udraws(:,t,n));
                    
                else
                    
                    for s = find(sNaN(:,t)')
                        ndxOther = 1 : Ns ~= s;
                        thisMu  = Sposterior(s) + beta1(s,:,t) * (S(ndxOther,t) - Sposterior(ndxOther));
                        thisSig = sqrtOmega1(s,t);
                        S(s,t)  = drawTruncNormal(thisMu, thisSig, elbBound, udraws(s,t,n));
                    end
                    
                end
            end
            
            
            Y(ndxS,t)   = S(:,t); % note: future values in YY need not be updated here, since we are looping forward
            Ytilde(:,t) = Y(:,t) - Y0(:,t);
            
        end % any(sNaN(:,t))
        
        % update (even if not(sNaN))
        if t >= p
            this            = Ytilde(:,t:-1:t-p+1);
            STATElag        = this(:);
        else
            this     = STATElag(1 : Ny * (p - 1));
            STATElag = cat(1, Ytilde(:,t), this);
        end
    end % for t
    
    if n > burnin
        shadowrateDraws(:,:,n-burnin) = S;
    end
    
    if showProgress
        progressbar(n / (totalNdraws))
    end
end


