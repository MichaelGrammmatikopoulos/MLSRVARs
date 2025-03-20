function diagonal = getmaindiagonal(A)
%GETMAINDIAGONAL: get the main diagonal of a square matrix 
% https://www.mathworks.com/matlabcentral/answers/355538-how-to-extract-the-diagonal-of-a-given-matrix

if nargin<1
    error('Too few inputs');
end

if nargin>1
    error('Too many inputs');
end

[i,j] = size(A);

if i~=j 
    error('matrix is not square');
else
    if i == 1 
        error('input is not a matrix');
    end
end

diagonal = A(sub2ind(size(A),1:size(A,1),1:size(A,2)));

end