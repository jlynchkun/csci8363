%Use matrix factorization to learn P in such a way that the rows of P 
%correspond to subjects and the columns of p correspond to k "discovered" 
%latent explanatory vectors
%
%Then for each discovered column, rank the values in the column, select 
%the top and bottom 50 or 100, and then do test to see if the group's 
%survival is statistically different.%

%B is the key - this is the gene/image correlation matrix

%%
%%load data
load Breast_cancer_image  %breast cancer image and gene data
load ppiMatrixTF  %gene names and gene associations database
load Breast_ge_image_corr_data_0_5 %correlation B matrix between genes and images - Pearson's rho
load image_name %names of image features - NEEDED?

%if gene expression data isn't in the interaction database, leave it out
Gene_expression_data = Gene_expression_data(:,~cellfun('isempty',GeneName(:,2))); 
Selected_GeneName = GeneName(~cellfun('isempty',GeneName(:,2)),2);
[Common_gene ai bi] = intersect(Selected_GeneName, gene_name_ge);
Gene_expression_data = Gene_expression_data(:,ai);

%sparse gene interaction network of labeled genes
gene_network = sparse(ppiMatrixTF(bi,bi));

%absolute value of correlations between image data and labeled genes
corr_network = abs(corr_data_0_5(ai,:));% - THIS IS B


%%
%INITIALIZE VARIABLES - make all gene and image data positive
disp('Initializing variables...')
X1 = Gene_expression_data;
X2 = Image_data;
X1 = [max(X1,0) max(-X1,0)];
X2 = [max(X2,0) max(-X2,0)];


%this is gene-gene association matrix (4 of these b/c we're doubling the X data to make positive/negative different)
A =[gene_network gene_network;
   gene_network gene_network];

%this is B
B =[corr_network corr_network;
    corr_network corr_network];

%% %Set parameters for Non-negative matrix factorization
L1 = 1; L2 = 1; %parameter for H1*A*H1^T - corresponds to LAMBDA_1
r1 = 1; %limit growth of W - corresponds to GAMMA_1
r2 = 1;%limit growth of H1 and H2
r3= 10;  %not used for anything as far as I can tell
K = 100;  %%%%K is the number of dimensions we're going to look at in our factorization

[n m1] = size(X1);
[n m2] = size(X2);



num_iterations = 100;
results = zeros(num_iterations*10,6);
results_counter = 1;
for c_cutoff = 0:.1:.9
corr_network(corr_network < c_cutoff) = 0; %this only works b/c i'm going increasing
results_temp = zeros(num_iterations,6);
    parfor iter=1:num_iterations

    
    %% initialize random factors for |X - WH1| + |X - WH2|
    %KEY - we only want to 
    W0 = rand(n,K);
    Initial_W = W0;  %W0 is initialized to random numbers between 0 and 1

    H10 = rand(K,m1);   %initialize H to random numbers between 0 and 1
    H20 = rand(K,m2);


    %% Call TaeHyun's functions for NNMF to get factorizations 
    %1 - without B
    %2 - with B
    %3 - LEARN B (introducing the 
    opt = 5; 

    [Initial_W,W1,H11,H12,W2,H21,H22,W3,H31,H32,W_ge,H_ge,W_image,H_image] = DNMF_residue_comodule_image_ge_assoc(X1,X2,W0,H10,H20, A,abs(B),L1,L2,r1,r2,r3, K, opt);



    %% Kaplan-Meyer Survival Analysis
    num_p = 50; %Number of patients to test in each group

    % CASE 1
    W = W1; %W from NOT B 
    pval = [];
    for i=1:K
        [val idx] = sort(W(:,i),'descend');[junka junkb junkc] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
        pval = [pval; junkc];
    end


    %CASE 2
    % Kaplan Meir survival analysis
    W = W2;
    pval2 = [];
    for i=1:K
        [val idx] = sort(W(:,i),'descend');[junka junkb junkc] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
        pval2 = [pval2; junkc];
    end


    %CASE 3
    % Kaplan Meir survival analysis
    W = W3;
    pval3 = [];
    for i=1:K
        [val idx] = sort(W(:,i),'descend');[junka junkb junkc] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
        pval3 = [pval3; junkc];
    end


    %CASE 4 - GENES ONLY
    W = W_ge;
    pval_ge = [];
    for i=1:K
        [val idx] = sort(W(:,i),'descend');[junka junkb junkc] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
        pval_ge = [pval_ge; junkc];
    end


    %CASE 5 - Image data only
    W = W_image;
    pval_image = [];
    for i=1:K
        [val idx] = sort(W(:,i),'descend');[junka junkb junkc] = logrank_no_fig(survival(idx(1:num_p)), survival(idx(end-num_p+1:end)));
        pval_image = [pval_image; junkc];
    end

    
    %RECORD RESULTS
    results_temp(iter,:) =     [c_cutoff length(find(pval < 0.05))  length(find(pval2 < 0.05)) length(find(pval3 < 0.05)) length(find(pval_ge < 0.05)) length(find(pval_image < 0.05))];

    end
    %combine results
    for ii=1:num_iterations
        results(results_counter,:) = results_temp(ii,:);
        results_counter = results_counter+1;
    end
end



