function [W, H] = nnmf_by_taehyun(X, W, H, iterations)
%NNMF - non-negative matrix factorization
%
% [W, H] = nnmf(V, r, iterations)
%
% Input:
%   V   - the matrix to factorize
%   r   - number of basis vectors to generate
%   iterations - number of EM iterations to perform
%
% Results:
%   W   - a set of r basis vectors
%   H   - represenations of the columns of V in 
%         the basis given by W
%
%   Author: David Ross
%   Follows: D.D. Lee & H.S. Seung, "Learning the parts of 
%            objets by non-negative matrix factorization"
%            NATURE, vol 401, 21 Oct. 1999.

% %---------------------------------------
% % check the input
% %---------------------------------------
% error(nargchk(3,3,nargin));
% 
% if ndims(V) ~= 2
%     error('V must have exactly 2 dimensions');
% end
% 
% if prod(size(r)) ~= 1 | r < 1
%     error('r must be a positive scalar');
% end
% 
% if prod(size(iterations)) ~= 1 | iterations < 1 
%     error('number of iterations must be a positive scalar');
% end
% 
%---------------------------------------
% Initialization
%---------------------------------------
V = X;
N = size(V,1); % dimensionality of examples (# rows)
C = size(V,2); % number of examples (columns)
% 
% W = 2 * rand(N,r);
% H = 2 * rand(r,C);

V = V + eps;
err = [];

%---------------------------------------
% EM
%---------------------------------------

for it_number = 1:iterations
   
    %% E-Step
    V_over_WH = V ./ (W * H + (V==0));
    W = W .* (V_over_WH * H');
    W = W ./ repmat(sum(W,1), [N 1]);
    
    %% M-Step
    V_over_WH = V ./ (W * H + (V==0));
    H = H .* (W' * V_over_WH);
    
    err = [err norm(V - W*H, 2)];
    
end