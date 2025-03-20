function ydraw = drawTruncNormal(mu,sqrtVCV,elb,rndStream)
% generates 1 draw from truncated scalar normal N(mu, sqrtVCV * sqrtVCV')
% subject to ydraws <= elb
% returns scalar ydraws 
%


%   Coded by  Elmar Mertens, em@elmarmertens.com

if ~isscalar(mu) && ~isscalar(sqrtVCV)
    error('this function is for scalars only')
end

if nargin < 4
    rndStream = getDefaultStream;
end
if isnumeric(rndStream)
    udraw   = rndStream;
else
    udraw   = rand(rndStream);
end

tol = 1e-10; % to decide when sqrtVCV diagonals are zero

N      = 1; 
if N ~= length(mu)
    error('dimension mismatch')
end


% make sure that sqrtVCV has positive diagonals
% (sqrt factors obtained from QR can have negative diagonals, for example)
sqrtVCV = abs(sqrtVCV);


thismu       = mu;

if sqrtVCV > tol
    ub      = (elb - thismu) ./ sqrtVCV;
    
    % PHIbar  = normcdf(ub);
    PHIbar  = 0.5 * erfc(- sqrt(0.5) * ub);
    
    % zdraw = norminv(udraw * PHIbar);
    if PHIbar > eps
        zdraw = -sqrt(2) * erfcinv(2 * udraw * PHIbar);
    else
        zdraw = ub;
        %     jack  = trandn(-Inf, ub);
    end

    ydraw  = thismu + sqrtVCV * zdraw;
    
else
    ydraw  = thismu;
end
