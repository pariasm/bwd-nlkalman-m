function d = L1_distance(a,b);
% L1_DISTANCE - computes Euclidean distance matrix
%
% E = L1_distance(A,B)
%
%    A - (DxM) matrix 
%    B - (DxN) matrix
%    df = 1, force diagonals to be zero; 0 (default), do not force
% Returns:
%    E - (MxN) L1 distances between the columns of A and B

% if (nargin < 3)
%    df = 0;    % by default, do not force 0 on the diagonal
% end

% if (size(a,1) ~= size(b,1))
%    error('A and B should be of same dimensionality');
% end

% if ~(isreal(a)*isreal(b))
%    disp('Warning: running distance.m with imaginary numbers.  Results may be off.'); 
% end
% 
if (size(a,1) == 1)
   a = [a; zeros(1,size(a,2))]; 
   b = [b; zeros(1,size(b,2))]; 
end

d=abs(repmat(a,[1 size(b,2)]) - b);
d = sum(d); 

% force 0 on the diagonal? 
% if (df==1)
%   d = d.*(1-eye(size(d)));
% end
