function [Xdraw, CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx] = ...
    VARTVPSVprecisionsamplerNaN(aaa,invbbb,y,yNaN,y0,ybar,CC,QQ,RR1,...
    arows, acols, asortndx, brows, bcols, bsortndx, Ndraws)
% VARTVPSVprecisionsamplerNaN0const for case of VAR(p) with missing values and fixed initial conditions
%
% USAGE: [Xdraw, CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx] = ...
%     VARTVPSVprecisionsamplerNaN(aaa,invbbb,ccc,y,yNaN,y0,ybar,rndStream,CC,QQ,RR1,...
%     arows, acols, asortndx, brows, bcols, bsortndx)
%
% aaa is Ny x Ny x p (x T) (T dimension is optional)
% invbbb is Ny x Ny (x T) (T dimension is optional); invbbb * invbbb' is inverse of variance of VAR residuals
% ccc is measurement vector (typically be identity matrix), in general Ny x Ny ( x T)
%
% y is Ny x T data matrix, and yNaN is Ny x T logical matrix that designates the missing values 
% values of y(yNaN) will be ignored
%
% y0 is Ny x p matrix of initial conditions (p lags of y)
% ybar is a vector of intercepts (Ny x 1) or a matrix of time-varying intercepts
% arguments after rndStream can be empty and will be returned as outputs for use in future calls
% Xdraws is Ny * T vector output (can be shaped to Ny x T)

%% VERSION INFO
% AUTHOR    : Elmar Mertens

if nargin < 16
    Ndraws = 1;
end

% get dimensions
[Ny, T] = size(y);
p       = size(aaa,3);
Ny0     = size(y0,1);
Nw      = size(invbbb,2);

if Nw ~= Ny
    error('Expecting Nw identical to Ny')
end
if Ny ~= Ny0
    error('Expecting Ny identical to Ny0')
end

if nargin < 8
    CC  = [];
    QQ  = [];
    RR1 = [];
    [arows, acols, asortndx, brows, bcols, bsortndx] = deal([]);
end

if ndims(aaa) <= 3
    aaa = repmat(aaa, [1 1 1 T]);
end
if ismatrix(invbbb)
    invbbb = repmat(invbbb, [1 1 T]);
end

Nx    = Ny; % for ease of comparison against more general state space routines

NyT   = Ny * T;
NwT   = Nw * T;
NxT   = Nx * T;


%% construct vectorized state space
Y     = reshape(y, NyT, 1);
Ynan  = reshape(yNaN, NyT, 1);
Y     = Y(~Ynan);

%% vectorize input matrices
NxNx         = Nx * Nx;
NxNxT        = NxNx * T;
invbbb       = reshape(invbbb, NxNxT, 1);

if isvector(ybar) && (length(ybar) ==  Ny)
    XX0  = repmat(ybar, T, 1);
else
    XX0  = ybar(:); 
end

%% adjust XX0 for initial conditions
% adjust for initial conditions
for k = 1 : min(p,T)
    xndx                       = (k-1) * Ny + (1 : Ny);
    theseInitialLags           = p - k + 1;
    thisA                      = reshape(aaa(:,:,k:p), Ny, Ny * theseInitialLags);
    thisX0                     = reshape(y0(:,1:p-k+1), Ny * theseInitialLags, 1);
    XX0(xndx)                  = XX0(xndx) + thisA * thisX0; 
end

%% CC and prepare Arows and Brows

if isempty(CC)

    % no pre-allocation of memory here, since to be evaluated only once

    % AA, build sequentially: first unit diagonal
    arows     = 1 : NxT;
    acols     = 1 : NxT;
    % add p lags (sequentially)
    for k = 1 : p
        theserows   = repmat((1 : Nx)', 1 , Nx, T - k);
        theserows   = theserows + permute(Nx * (k : T-1), [1 3 2]);
        arows       = [arows(:); theserows(:)];

        thesecols   = repmat(1 : Nx * (T - k), Nx, 1);
        acols       = [acols(:); thesecols(:)];
    end

    % sort A indices
    ndx = sub2ind([NxT, NxT], arows, acols);
    [~, asortndx] = sort(ndx);
    arows         = arows(asortndx);
    acols         = acols(asortndx);

    % BB
    brows  = repmat((1 : Nx)', 1 , Nw, T);
    brows  = brows + permute(Nw * (0 : T-1), [1 3 2]);
    brows  = reshape(brows, Nx * NwT, 1);
    
    bcols  = repmat(1 : NwT, Nx, 1);
    bcols  = reshape(bcols, NwT * Nx, 1);

    % sort B indices
    ndx = sub2ind([NxT, NwT], brows, bcols);
    [~, bsortndx] = sort(ndx);
    brows         = brows(bsortndx);
    bcols         = bcols(bsortndx);

    % C
    CC        = speye(NyT);
    % drop rows associated with NaN
    CC        = CC(~yNaN,:);
    % perform QR
    [QQ,RR]   = qr(CC');
    [N1, N2]  = size(CC);
    N2        = N2 - N1;
    RR1       = RR(1:N1,1:N1)';

else

    N1        = size(RR1,1);
    N2        = size(QQ,1) - N1;
    
end

QQ1       = QQ(:,1:N1)';
QQ2       = QQ(:,N1+1:end)';

%% sparse builds for BB and AA
% AA
pEff             = min(T,p); % to catch cases where T<p
Naa              = NxT + NxNx * p * max(T - p, 0) + sum(NxNx * (1 : pEff - 1));
values           = ones(Naa,1);
% fill p lags (sequentially)
offset = NxT;
for k = 1 : pEff
    values(offset + (1 : NxNx * (T-k))) = -reshape(aaa(:,:,k,1+k:T), Nx * Nx * (T - k), 1);
    offset = offset + NxNx * (T-k);
end
values              = values(asortndx);
AA                  = sparse(arows, acols, values, NxT, NxT);

% BB
invbbb  = invbbb(bsortndx);
invBB   = sparse(brows, bcols, invbbb, NxT, NxT);

%% means and innovations
EX        = AA \ XX0;
EY        = CC * EX;

X1tilde   = RR1 \ (Y - EY);

QQX1tilde = QQ1' * X1tilde;

%% precision-based sampler
AAtilde       = invBB * AA;
AAtildeQQX1   = AAtilde * QQX1tilde;
AAtildeQQ2    = AAtilde * QQ2';
invQSIG22     = transpose(AAtildeQQ2) * AAtildeQQ2;
cholinvQSIG22 = chol(invQSIG22, 'lower');

X2hat         = - cholinvQSIG22 \ (AAtildeQQ2' * AAtildeQQX1);

Z2draw        = randn(N2, Ndraws) + X2hat;
X2draw        = cholinvQSIG22' \ Z2draw;
Xdraw         = EX + QQX1tilde + QQ2' * X2draw;

