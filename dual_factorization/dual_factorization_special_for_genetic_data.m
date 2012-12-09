function [results] = dual_factorization(GeneName,Gene_expression_data,Image_data,survival, gene_name_ge,ppiMatrixTF,K)
%% Use matrix factorization to learn W in such a way that the rows of W 
%correspond to subjects and the columns of W correspond to k "discovered" 
%latent explanatory vectors
%Then for each discovered column, rank the values in the column, select 
%the top and bottom 50 or 100, and then do test to see if the groups' 
%survival rates are statistically different.%

%INPUT - X1, X2, labels, A matrix of interactions between X1 elements, k
%(number of dimensions for factorization)
%OUTPUT - a 5xk matrix of p-values by columns for each method showing
%significance for group distinctions using kaplan-meyer

%B is the key - this is the gene/image correlation matrix


%% Clean up gene data
[Common_gene ai bi] = intersect(GeneName(:,2), gene_name_ge);
Gene_expression_data = Gene_expression_data(:,ai);
preprocess = 0; %remove data with small absolute values and small variance

%%
%INITIALIZE VARIABLES - make all gene and image data positive
disp('Initializing variables...')
gene_network = sparse(ppiMatrixTF(bi,bi)); %genes from interaction network also labeled in breast cancer data
%NOTE: changed X1 and X2 to be consistent with Zhang paper
X1 = Image_data;

%preprocess genetic data - OPTIONAL - THIS WAS DONE IN ZHANG 2011
if preprocess
    mask1 = genelowvalfilter(Gene_expression_data','Percentile',60);
    mask2 = genevarfilter(Gene_expression_data','Percentile',30);    
    X2 = Gene_expression_data(:,mask1 & mask2); %I THINK this applies both masks
    gene_network = gene_network(mask1&mask2,mask1&mask2);
else
    X2 = Gene_expression_data;
end

%Create B - correlation matrix between Gene_expression_data and Image_data
%THIS SHOULD BE SIMILAR TO MULTIPLYING H1'*H2 - CHECK TO SEE IF IT IS
%QUESTION - SHOULD WE JUST DO THIS LIKE Zhang DOES?
disp('Calculating gene-image correlation matrix...')
correlations = corrcoef([X1,X2]);
B = correlations(end-size(X1,2)+1:end,1:size(X2,2)); %B is just the lower quarter of this correlation matrix - WEIRD DISTRIBUTION?
B = abs([B B; B B]); %double to make the same shape as new X's
B(B<.7)=0; %make sparse - DO THIS???
size(B)
%make all positive
X1 = [max(X1,0) max(-X1,0)]; %IMAGE DATA
X2 = [max(X2,0) max(-X2,0)]; %GENETIC DATA

%this is gene-gene association matrix (4 of these b/c we're doubling the X data to make positive/negative different)
A =[gene_network gene_network;
   gene_network gene_network];



[n,m1] = size(X1);
[n,m2] = size(X2);

%keep track of best factorizations for finding comodules
bestW_1=zeros(n,K);
bestH1_1=zeros(K,m1);
bestH2_1=zeros(K,m2);
bestobj1_1=1000000000;
bestobj2_1=1000000000;
bestW_2=zeros(n,K);
bestH1_2=zeros(K,m1);
bestH2_2=zeros(K,m2);
bestobj1_2=1000000000;
bestobj2_2=1000000000;
bestW_3=zeros(n,K);
bestH1_3=zeros(K,m1);
bestH2_3=zeros(K,m2);
bestobj1_3=1000000000;
bestobj2_3=1000000000;



%% Do updates
gamma1=20;
gamma2=10;
lambda1=.0001;
lambda2=.01;
lambda3=.0001;

%FOR LOOP GOES HERE
all_results = {};
for iter=1:5

% initialize random factors for |X - WH1| + |X - WH2|
%KEY - we only want to use a common W/H for each run to be sure that we
%are comparing apples to apples as far as results go
W = rand(n,K);
H1 = rand(K,m1);   %initialize H to random numbers between 0 and 1
H2 = rand(K,m2);

%% UPDATES
%case 1 - omit B, so lambda2 = 0, lambda3 = 0
[W1_1,H1_1,H2_1] = multiplicative_update(X1,X2,W,H1,H2,A,B,lambda1,0,0,gamma1,gamma2,K);

%case 2 - keep B, default from Zhang paper, don't learn
%B so lambda3 = 0
[W1_2,H1_2,H2_2] = multiplicative_update(X1,X2,W,H1,H2,A,B,lambda1,lambda2,0,gamma1,gamma2,K);

%case 3 - learn B, so lambda3 != 0 
[W1_3,H1_3,H2_3] = multiplicative_update(X1,X2,W,H1,H2,A,B,lambda1,lambda2,lambda3,gamma1,gamma2,K);

%case 4 - only image data
maxiter = 50;
[W1_image, H1_image] = non_negative_matrix_factorization(X1, W, H1, maxiter);

%case 5 - only genetic data
[W1_gene, H1_gene] = non_negative_matrix_factorization(X2, W, H2, maxiter);


%% SIGNIFICANCE OF RESULT
%Do logrank test (Mantel-Cox test) statistics
%In each case, go through each column and sort the
%values. Take the highest and lowest num_p subjects
%and see if their survival rates are significantly
%different using the Mantel-Cox test.  Record the
%p-values from the test results.  

num_p = 50; %Number of patients to test in each group

% CASE 1
pval1 = [];
for i=1:K
    [val idx] = sort(W1_1(:,i),'descend');[junka junkb p] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
    pval1 = [pval1; p];
end


%CASE 2
pval2 = [];
for i=1:K
    [val idx] = sort(W1_2(:,i),'descend');[junka junkb p] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
    pval2 = [pval2; p];
end


%CASE 3
pval3 = [];
for i=1:K
    [val idx] = sort(W1_3(:,i),'descend');[junka junkb p] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
    pval3 = [pval3; p];
end


%CASE 4 - Image only
pval_image = [];
for i=1:K
    [val idx] = sort(W1_image(:,i),'descend');[junka junkb p] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
    pval_image = [pval_image; p];
end


%CASE 5 - Gene data only
pval_gene = [];
for i=1:K
    [val idx] = sort(W1_gene(:,i),'descend');[junka junkb p] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
    pval_gene = [pval_gene; p];
end
%RECORD RESULTS
%results are p-values for each method, with each column corresponding to a
%different method used.
results = [pval1; pval2; pval3; pval_image; pval_gene];


%% Comodules

    % compute residue
    %method 1
    newobj1_1 = sum(sum((X1-W1_1*H1_1).^2));
    newobj2_1 = sum(sum((X2-W1_1*H2_1).^2));    
    if (newobj1_1<bestobj1_1)||(newobj2_1<bestobj2_1)        
        bestobj1_1 = newobj1_1;
        bestobj2_1 = newobj2_1;  
        bestW1_1 = W1_1;
        bestH1_1 = H1_1;
        bestH2_1 = H2_1;
    end
    %method 2
    newobj1_2 = sum(sum((X1-W1_2*H1_2).^2));
    newobj2_2 = sum(sum((X2-W1_2*H2_2).^2));    
    if (newobj1_2<bestobj1_2)||(newobj2_2<bestobj2_2)        
        bestobj1_2 = newobj1_2;
        bestobj2_2 = newobj2_2;  
        bestW1_2 = W1_2;
        bestH1_2 = H1_2;
        bestH2_2 = H2_2;
    end
    %method 3
    newobj1_3 = sum(sum((X1-W1_3*H1_3).^2));
    newobj2_3 = sum(sum((X2-W1_3*H2_3).^2));    
    if (newobj1_3<bestobj1_3)||(newobj2_3<bestobj2_3)        
        bestobj1_3 = newobj1_3;
        bestobj2_3 = newobj2_3;  
        bestW1_3 = W1_3;
        bestH1_3 = H1_3;
        bestH2_3 = H2_3;
    end

all_results{iter} = results;

end
%END FOR LOOP for iterations


%% get comodules(index)
tt = 7;
isdouble=1;
[Co_module_1] = DNMF_comodule(bestW1_1, bestH1_1, bestH2_1, tt, isdouble);
[Co_module_2] = DNMF_comodule(bestW1_2, bestH1_2, bestH2_2, tt, isdouble);
[Co_module_3] = DNMF_comodule(bestW1_3, bestH1_3, bestH2_3, tt, isdouble);



save(['factors_5_iterations_high_threshold_k_' num2str(K) '.mat'],'Co_module_1','Co_module_2','Co_module_3','bestW1_1','bestW1_2','bestW1_3','W1_image','W1_gene','bestH1_1','bestH1_2','bestH1_3','H1_image','H1_gene','bestH2_1','bestH2_2','bestH2_3','all_results');