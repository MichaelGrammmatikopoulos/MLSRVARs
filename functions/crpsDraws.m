function crps = crpsDraws(y, Xdraws, hasBeenSorted)
% CRPSDRAWS ...
%
% implements Krueger at el efficient formula
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 15-Jun-2016 04:40:44 $
% $Revision : 1.00 $
% DEVELOPED : 8.6.0.267246 (R2015b)
% FILENAME  : crpsDraws.m

Xdraws    = Xdraws(:);
Ndraws    = length(Xdraws);

if any(isnan(Xdraws)) || isnan(y)
    crps  = NaN;
    return;
end

% nbins         = ceil(Ndraws / 2);

if nargin < 3
    hasBeenSorted = [];
end
if isempty(hasBeenSorted)
    hasBeenSorted = issorted(Xdraws);
end


%% sort if needed
if hasBeenSorted
    Xordered  = Xdraws;
else
    Xordered = sort(Xdraws);
end

%% compute CRPS

[alpha, beta] = deal(NaN(Ndraws-1,1)); % for cases i=1:end-1

% cases i=1:end-1
ndx1 = y < Xordered(1:end-1);
ndx2 = (Xordered(1:end-1) <= y) & (y <= Xordered(2:end));
ndx3 = Xordered(2:end) < y;

alpha(ndx1) = 0;

foo = y - Xordered(1:end-1);
alpha(ndx2) = foo(ndx2);

foo = Xordered(2:end) - Xordered(1:end-1);
alpha(ndx3) = foo(ndx3);

foo        = Xordered(2:end) - Xordered(1:end-1);
beta(ndx1) = foo(ndx1);

foo = Xordered(2:end) - y;
beta(ndx2) = foo(ndx2);
beta(ndx3) = 0;

% case i = 0 and end
ALPHA = zeros(Ndraws+1,1);
BETA  = ALPHA;

ALPHA(2:end-1) = alpha;
BETA(2:end-1) = beta;
if y >= Xordered(end)
    ALPHA([1 end]) = y - Xordered(end);
end
if y <= Xordered(1)
    BETA([1 end]) = Xordered(1)-y;
end

ndx  = (0:Ndraws)' / Ndraws;
crps = sum(ALPHA .* ndx.^2 + BETA .* (1 - ndx).^2);