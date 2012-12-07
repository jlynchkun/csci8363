% synthetic data
% need survival data
s_survival = [1;
              0]

% image data
s_X1 = [2 1 -2 0;
        1 1 0 -2] + 0.01*randn(2,4)

s_X1_z = zscore(s_X1)
      
s_X1_nn = [max(s_X1_z,0) max(-s_X1_z,0)] %IMAGE DATA

% gene expression data
% two subjects!
s_X2 = [1 -1 0;
        1 0 1] + 0.01*randn(2,3)

s_X2_z = zscore(s_X2)
s_X2_nn = [max(s_X2_z,0) max(-s_X2_z,0)] %GENETIC DATA


% gene-gene interaction
s_gene_network = [0 1 0;
                  0 0 1;
                  1 0 0]
     
s_A =[s_gene_network s_gene_network;
    s_gene_network s_gene_network];
  
correlations = corrcoef([s_X1_nn,s_X2_nn])
s_B = correlations(end-size(s_X1_nn,2)+1:end,1:size(s_X2_nn,2)) %B is just the lower quarter of this correlation matrix - WEIRD DISTRIBUTION?
s_B = abs(s_B)


