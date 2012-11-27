function [W,H1,H2] = DNMF_mm(X1, X2, A, B, r1, r2, L1, L2, K, maxiter, speak, fid)
%
% Multiple NMF using euclidean distance update equations:
%
% Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix
% Factorization', Adv. Neural Info. Proc. Syst. 13, 556-562.
%
% INPUT:
% X1 (N,M1): N (dimensionallity) x M1 (feature 1) non negative input matrix
% X2 (N,M2): N (dimensionallity) x M2 (feature 2) non negative input matrix
% A (M2,M2): M2 x M2, non negative input matrix
% B (M1,M2): M1 x M2, non negative input matrix 

% r1       : limit the growth of W
% r2       : limit the growth of H1 and H2
% L1       : weigh the must link constraints in A
% L2       : weigh the must link constraints in B
% K        : Number of components
% maxiter  : Maximum number of iterations to run
% speak    : prints iteration count and changes in connectivity matrix 
%            elements unless speak is 0
%fid       : file identifier which is used to store the changes record

% OUTPUT:
% W        : N x K matrix
% H1       : K x M1 matrix
% H2       : K x M2 matrix
%
% Shihua Zhang
% Computational Biology program in the Department of Biological Sciences
% University of Southern California
% zsh@amss.ac.cn
% 2009/05/10

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print_iter = 1; % iterations between print on screen and convergence test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for negative values in X1 and X2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (min(min(X1)) < 0) || (min(min(X2)) < 0) 
    error('Input matrix elements can not be negative');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for same rows in X1 and X2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n1,m1] = size(X1);
[n2,m2] = size(X2);

if n1 ~= n2 
    error('Input matrices should have the same rows');
    return
end
n = n1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize random W, H1 and H2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = rand(n,K);
H1 = rand(K,m1);
H2 = rand(K,m2);


% use W*H to test for convergence
Xr_old1 = W*H1;
Xr_old2 = W*H2;

for iter = 1 : maxiter
    % Euclidean multiplicative method
    W = W.*([H1 H2]*[X1 X2]')'./(W*([H1 H2]*[H1 H2]'+r1*eye(K))+eps);    %Update rule-1;   
  
    HH1 = H1.*(W'*X1 + L2/2*H2*B')./((W'*W+r2*ones(K))*H1+eps);     %%Update H1 H2 simultaneously
    H2 = H2.*(W'*X2 + L1*H2*A + L2/2*H1*B)./((W'*W+r2*ones(K))*H2+eps); 
    H1 = HH1;

   % iter
    % print to screen
    if (rem(iter,print_iter) == 0) & speak
        Xr1 = W*H1+eps;            
        Xr2 = W*H2+eps;
        diff_step = sum(sum(abs(Xr_old1-Xr1)))+sum(sum(abs(Xr_old2-Xr2)));
        
        diff2 = -L1*trace(H2*A*H2');
        diff3 = -L2*trace(H1*B*H2');
        diff4 = r1*sum(sum(W.*W));
	    diff5 = r2*(sum(sum(H1).^2)+sum(sum(H2).^2));
       
	    Xr_old1 = Xr1;
        Xr_old2 = Xr2;
        eucl_dist1 = nmf_euclidean_dist(X1,W*H1);
        eucl_dist2 = nmf_euclidean_dist(X2,W*H2);
        diff1 = eucl_dist1 + eucl_dist2;
		diff = diff1 + diff2 + diff3 + diff4 + diff5;
        errorx1 = mean(mean(abs(X1-W*H1)))/mean(mean(X1));
        errorx2 = mean(mean(abs(X2-W*H2)))/mean(mean(X2));
        errorx = errorx1 + errorx2;
        fprintf(fid,'%s\n',[sprintf('Iter = \t'),int2str(iter),...
            sprintf('\t relative error = \t'),num2str(errorx),...
            sprintf('\t diff_step = \t'),num2str(diff_step),...
            sprintf('\t diff = \t'), num2str(diff),...
            sprintf('\t diff1 = \t'), num2str(diff1),...
            sprintf('\t diff2 = \t'), num2str(diff2),...
            sprintf('\t diff3 = \t'), num2str(diff3),...
		    sprintf('\t diff4 = \t'), num2str(diff4),...
		    sprintf('\t diff5 = \t'), num2str(diff5)]);
        if errorx < 10^(-5), break, end
    end
end

function err = nmf_euclidean_dist(X,Y)
err = sum(sum((X-Y).^2));
