function [return1 return2] = main(GeneName,Gene_expression_data,Image_data,survival, gene_name_ge,ppiMatrixTF)
%% Use matrix factorization to learn W in such a way that the rows of W 
%correspond to subjects and the columns of W correspond to k "discovered" 
%latent explanatory vectors
%
%Then for each discovered column, rank the values in the column, select 
%the top and bottom 50 or 100, and then do test to see if the groups' 
%survival rates are statistically different.%

%to use, first 
%load Breast_cancer_image
%load ppiMatrixTF

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


%make all positive
X1 = [max(X1,0) max(-X1,0)]; %IMAGE DATA
X2 = [max(X2,0) max(-X2,0)]; %GENETIC DATA

%this is gene-gene association matrix (4 of these b/c we're doubling the X data to make positive/negative different)
A =[gene_network gene_network;
   gene_network gene_network];

%Create B - correlation matrix between Gene_expression_data and Image_data
%THIS SHOULD BE SIMILAR TO MULTIPLYING H1'*H2 - CHECK TO SEE IF IT IS
%QUESTION - SHOULD WE JUST DO THIS LIKE Zhang DOES?
disp('Calculating gene-image correlation matrix...')
correlations = corrcoef([X1,X2]);
B = correlations(end-size(X1,2)+1:end,1:size(X2,2)); %B is just the lower quarter of this correlation matrix - WEIRD DISTRIBUTION?
B = abs(B);






%% Do updates
num_iterations = 1;
results_counter = 1;

%dimensionality of result
K=50;

%range of parameters to be tested
parameter_values = [.0001 .001 .01 .1 1 10 20];

%results are [B_cutoff lambda_1 lambda_2 gamma_1 gamma_2 case1 case2 case3 case4 case5]
%where "case#" is the number of significant columns found in W with regard
%to the Kaplam-meyer statistic
results = [];

for gamma1 = parameter_values
    for gamma2 = parameter_values
        for lambda1 = parameter_values
            for lambda2 = parameter_values
                for lambda3 = parameter_values
                    gamma1
                    gamma2
                    lambda1
                    lambda2
                    lambda3
                    results_temp = [];
                    for iter=1:num_iterations
                        iter
                        % initialize random factors for |X - WH1| + |X - WH2|
                        %KEY - we only want to use a common W/H for each run to be sure that we
                        %are comparing apples to apples as far as results go
                        [n,m1] = size(X1);
                        [n,m2] = size(X2);
                        W = rand(n,K);
                        H1 = rand(K,m1);   %initialize H to random numbers between 0 and 1
                        H2 = rand(K,m2);

                        %% UPDATES
                        %case 1 - omit B, so gamma2 = 0, lambda3 = 0
                        [W1_1,H1_1,H2_1] = multiplicative_update(X1,X2,W,H1,H2,A,B,lambda1,lambda2,0,gamma1,0,K);

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
                        sig_cutoff = 0.05;
                        results_temp(iter,:) = [gamma1 gamma2 lambda1 lambda2 lambda3 length(find(pval1 < sig_cutoff))  length(find(pval2 < sig_cutoff)) length(find(pval3 < sig_cutoff)) length(find(pval_image < sig_cutoff)) length(find(pval_gene < sig_cutoff))];

                        

                        
                    end
                    %combine results
                    for ii=1:num_iterations
                        results(results_counter,:) = results_temp(ii,:);
                        results_counter = results_counter+1;
                    end   
                    
                    
                    
                end
            end
        end
    end
end

% 
% %% %Set parameters for Non-negative matrix factorization
% L1 = 1; L2 = 1; %parameter for H1*A*H1^T - corresponds to LAMBDA_1
% r1 = 1; %limit growth of W - corresponds to GAMMA_1
% r2 = 1;%limit growth of H1 and H2
% K = 100;  %%%%K is the number of dimensions we're going to look at in our factorization
% 
% [n m1] = size(X1);
% [n m2] = size(X2);
% 
% 
% 
% 
% for c_cutoff = 0:.1:.9
% corr_network(corr_network < c_cutoff) = 0; %this only works b/c i'm going increasing
% results_temp = zeros(num_iterations,6);
%     parfor iter=1:num_iterations
% 
%     
%     %% initialize random factors for |X - WH1| + |X - WH2|
%     %KEY - we only want to 
%     W0 = rand(n,K);
%     Initial_W = W0;  %W0 is initialized to random numbers between 0 and 1
% 
%     H10 = rand(K,m1);   %initialize H to random numbers between 0 and 1
%     H20 = rand(K,m2);
% 
% 
%     %% Call TaeHyun's functions for NNMF to get factorizations 
%     %1 - without B
%     %2 - with B
%     %3 - LEARN B (introducing the B-B0 term) 
%     %4 - Gene data only
%     %5 - Image data only
% 
%     [Initial_W,W1,H11,H12,W2,H21,H22,W3,H31,H32,W_ge,H_ge,W_image,H_image] = DNMF_residue_comodule_image_ge_assoc(X1,X2,W0,H10,H20, A,abs(B),L1,L2,r1,r2,r3, K, opt);
% 
% 
% 
%     %% Kaplan-Meyer Survival Analysis
%     num_p = 50; %Number of patients to test in each group
% 
%     % CASE 1
%     W = W1; %W from NOT B 
%     pval = [];
%     for i=1:K
%         [val idx] = sort(W(:,i),'descend');[junka junkb junkc] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
%         pval = [pval; junkc];
%     end
% 
% 
%     %CASE 2
%     % Kaplan Meir survival analysis
%     W = W2;
%     pval2 = [];
%     for i=1:K
%         [val idx] = sort(W(:,i),'descend');[junka junkb junkc] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
%         pval2 = [pval2; junkc];
%     end
% 
% 
%     %CASE 3
%     % Kaplan Meir survival analysis
%     W = W3;
%     pval3 = [];
%     for i=1:K
%         [val idx] = sort(W(:,i),'descend');[junka junkb junkc] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
%         pval3 = [pval3; junkc];
%     end
% 
% 
%     %CASE 4 - GENES ONLY
%     W = W_ge;
%     pval_ge = [];
%     for i=1:K
%         [val idx] = sort(W(:,i),'descend');[junka junkb junkc] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
%         pval_ge = [pval_ge; junkc];
%     end
% 
% 
%     %CASE 5 - Image data only
%     W = W_image;
%     pval_image = [];
%     for i=1:K
%         [val idx] = sort(W(:,i),'descend');[junka junkb junkc] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
%         pval_image = [pval_image; junkc];
%     end
% 
%     
%     %RECORD RESULTS
%     results_temp(iter,:) =     [c_cutoff length(find(pval < 0.05))  length(find(pval2 < 0.05)) length(find(pval3 < 0.05)) length(find(pval_ge < 0.05)) length(find(pval_image < 0.05))];
% 
%     end
%     %combine results
%     for ii=1:num_iterations
%         results(results_counter,:) = results_temp(ii,:);
%         results_counter = results_counter+1;
%     end
% end



