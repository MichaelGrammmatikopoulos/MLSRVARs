function r1 = rowextract(Xmat,which_row)
%COLUMNEXTRACT Summary of this function goes here
%   Detailed explanation goes here
if nargin<1
    error('too few inputs');
end
if nargin > 2
    error('too many arguments (need only matrix and row id)');
end
if ~isa(which_row,'double')
    error('second input is not a double: you need to identify which row to extract');
end
r1 = Xmat(which_row,:);
end