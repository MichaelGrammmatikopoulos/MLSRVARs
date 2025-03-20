function c1 = columnextract(Xmat,which_column)
%COLUMNEXTRACT Summary of this function goes here
%   Detailed explanation goes here
if nargin<1
    error('too few inputs');
end
if nargin > 2
    error('too many arguments (need only matrix and column id)');
end
if ~isa(which_column,'double')
    error('second input is not a double: you need to identify which column to extract');
end

c1 = Xmat(:,which_column);

end