function [G1, G2, P,error] = gene_image_optimization(Gene_expression_data,Image_data,labels,num_iterations)

%X1 = Gene expression data (n x m)
%X2 = Image Features (n x l)
%P = discovered patient groups (n x 2)
%P0 = initial patient groups (metastasis, survival, etc)
%G1 = meta gene expression (to be discovered) <- weights to patient groups
%   to create gene features X1
%G2 = meta image features (to be discovered) <- weights to patient groups
%   to create image features X2

%goal - to do non-negative matrix factorizations of X1 and X2 and minimize
%|X1 - P*G1|_F + |X2 - P*G2|_F + |P - P0|
%n x m  -  (n x 2)*(2 x m)


%Approach: Fix P (initialize as P0), find G1 and G2. Then fix G1 and G2, 
%find P.  Continue until approximation converges.

%norm(X,'fro') is the frobenius norm



%initialize variables
%put everything in the right order
[Y I] = sort(labels,'descend'); %sort labels by class - I is indices
labels_sorted = labels(I);
X1 = Gene_expression_data(I,1:3000);
X2 = Image_data(I,1:3000);
P0 = [labels_sorted > 0 labels_sorted < 0];


%make X1 and X2 into [pos neg] matrices, so there are all positive values
X1 = [(X1 .* (X1 > 0))   (-1*X1 .* (X1<0))];
X2 = [(X2 .* (X2 > 0))   (-1*X2 .* (X2<0))];


%sizes
size_n = size(X1,1) %number of patients
size_m = size(X1,2) %number of gene values
size_l = size(X2,2) %number of image features



%% experiment - see if sparsity makes this run better.
%UNFORTUNATELY IT DOESN'T SEEM TO - still too big for memroy, even with 9(%
%sparsity. However, just using a subset of the data columns DOES work, so
%perhaps find some better subsets (ones correlated with top singular
%values??  3000 columns each seems to work fi
% 
% %sparsity proportion
% s = 0.99;
% nnz(X1)
% sparse_X1 = rand(size(X1))>s;
% sparse_X2 = rand(size(X2))>s;
% X1 = X1 & sparse_X1;
% X2 = X2 & sparse_X2;
% nnz(X1)




%NOTE:
%matlab's NNMF isn't great for this because there are TWO matrices that
%need to be factored, but one of the factors is the same (namely P)

%fix P, solve for G. Then Fix G, solve for P. Repeat.
%Is this better than solving all at once? Yes - CVX doesn't allow 
%multiplication between teo variables that are being minimized

error = zeros(1,num_iterations);

%fix P initially
P = P0;
for i=1:num_iterations
i
    %P fixed, solve for G
    cvx_begin
        %declare variables
        variable G1( [2 size_m] );
        variable G2( [2 size_l] );

        minimize( norm(X1 - P*G1,'fro') + norm(X2 - P*G2,'fro') + norm(P-P0,2) );

        subject to
            G1 >= 0;
            G2 >= 0;

    cvx_end
    
    %fix G, solve for P
    cvx_begin
        variable P( [size_n 2])
        
        minimize( norm(X1-P*G1,'fro') + norm(X2 - P*G2,'fro') + norm(P-P0,2) )
        
        %subject to
        %    P >= 0 %is this necessary?
            
    cvx_end

    %record error
    error(i) = norm(X1-P*G1,'fro') + norm(X2 - P*G2,'fro') + norm(P-P0,2);
end

              
end

        
        
        